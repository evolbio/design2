#include <cmath>
#include "Individual.h"
#include "Performance.h"

// declare static private member variables

double Individual::mut;
double Individual::rec;
int Individual::totalLoci;
Allele Individual::mutStep;
double Individual::aSD;
double Individual::fitVar;
double Individual::gamma;
Loop Individual::loop;
int Individual::mutLocus;
double Individual::stochWt;
bool Individual::stoch;
ulong Individual::negLog2Rec;

// Algorithm for fast Poisson for lambda < 30
// from https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/
// Test measurement suggests about twice as fast as rnd.poisson()

int MyRandomPoisson(double mean)
{
    int k = 0;
    double p = 1.0;
    double L = exp(-mean);
    while (p > L){
        k++;
        p *= rnd.rU01();
    }
    return k-1;
}

// copy constructor
Individual::Individual(const Individual& other)
{
    genotype = std::unique_ptr<Allele[]> {new Allele[totalLoci]};
    if (stoch) stochast = std::unique_ptr<Allele[]> {new Allele[totalLoci]};
    for (int i = 0; i < totalLoci; ++i){
        genotype[i] = other.genotype[i];
        if (stoch) stochast[i] = other.stochast[i];
    }
    fitness = calcFitness();
}

// assignment constructor
Individual& Individual::operator=(const Individual& other)
{
    genotype = std::unique_ptr<Allele[]> {new Allele[totalLoci]};
    if (stoch) stochast = std::unique_ptr<Allele[]> {new Allele[totalLoci]};
    for (int i = 0; i < totalLoci; ++i){
        genotype[i] = other.genotype[i];
        if (stoch) stochast[i] = other.stochast[i];
    }
    fitness = calcFitness();
    return *this;
}

void Individual::initialize()
{
    genotype = std::unique_ptr<Allele[]> {new Allele[totalLoci]};
    // init to zero with {} initializer
    if (stoch) stochast = std::unique_ptr<Allele[]> {new Allele[totalLoci]{}};
    float p = static_cast<float>(1.0/sqrt(gamma));
    // p0 = 0 by assumption
    genotype[2] = p;                // q0
    genotype[3] = static_cast<float>(sqrt(1+gamma)*p);  // q1
    genotype[4] = p;                // q2
    switch(loop){                   // p2
        case Loop::open:
            genotype[0] = 1.0;      // p1
            genotype[1] = p;        // p2
            break;
        case Loop::close:
            genotype[0] = 1.0;      // p1
            genotype[1] = 0.0;      // p2
            break;
        case Loop::dclose:
            Allele r = 10.0f;       // by assumption, see plasticity II
            Allele k = 1/r;         // by assumption, see plasticity II
            genotype[0] = k;        // p1
            genotype[1] = k*r-p;    // p2
            genotype[5] = r;        // r
            genotype[6] = k;        // k
            break;
    }

    if (mutLocus >= 0){
        genotype[mutLocus] = mutateStep(genotype[mutLocus]);
        if (stoch) stochast [mutLocus] = mutateStep(stochast[mutLocus]);
    }
    else{
        for (int i = 0; i < totalLoci; ++i){
            genotype[i] = mutateStep(genotype[i]);
            if (stoch){
                stochast[i] = mutateStep(stochast[i]);
                if (stochast[i] < 0) stochast[i] = static_cast<Allele>(0);
            }
        }
    }
    fitness = calcFitness();
}

// set static variables used by class

void Individual::setParam(Param& param)
{
    mut = param.mutation;
    rec = param.recombination;
    totalLoci = param.loci;
    mutStep = param.mutStep;
    aSD = param.aSD;
    fitVar = param.fitVar;
    gamma = param.gamma;
    loop = param.loop;
    mutLocus = param.mutLocus;
    stochWt = param.stochWt;
    stoch = param.stoch;
    negLog2Rec = 1;         // set elsewhere when needed, here is just default value
}

// could use bit cache for random bits to speed up

Allele Individual::mutateStep(Allele a)
{
    auto c = static_cast<Allele>(rnd.rUniform(-mutStep,mutStep));
    return a + c;
}

// mut is per genotype mutation rate

void Individual::mutate()
{
    mutateG(genotype, false);
    if (stoch) mutateG(stochast, true);
}

// s == true => set negative values to zero, used for stochastic parameters which are standard deviations and so must be nonnegative
void Individual::mutateG(std::unique_ptr<Allele []>& g, bool s)
{
    if (mutLocus >= 0){
        if (rnd.rU01() < mut)
            g[mutLocus] = mutateStep(g[mutLocus]);
    }
    else{
        int hits = MyRandomPoisson(mut*totalLoci);  // about twice as fast as rnd.poisson()
        for (int i = 0; i < hits; ++i){
            ulong locus = rnd.rtop(totalLoci);
            g[locus] = mutateStep(g[locus]);
            if (s && (g[locus] < 0)) g[locus] = static_cast<Allele>(0);
        }
    }
}

// Possible modification: For example, if rec = 0.4, then use 10/25 to test, ie, use random numbers on [0..24] and success if random number is <= 9. Not sure how much performance boost would be obtained. Would only need five bits to generate random number, so could use bit operations on 64 bit random number.

void SetBabyGenotype(Individual& Parent1, Individual& Parent2, Individual& baby)
{
    auto& g1 = Parent1.genotype;
    auto& g2 = Parent2.genotype;
    auto& gb = baby.genotype;
    auto& s1 = Parent1.stochast;
    auto& s2 = Parent2.stochast;
    auto& sb = baby.stochast;
    ulong chrFlag = rnd.rbit();        // determines which parent is used for copying

    for (int i = 0; i < Parent1.totalLoci; ++i){
        gb[i] = (chrFlag) ? g1[i] : g2[i];
        if (Individual::stoch) sb[i] = (chrFlag) ? s1[i] : s2[i];
        if (rnd.rU01() < baby.rec) chrFlag ^= 1;        // flip flag if recombination at rate 0.5
    }
    baby.fitness = baby.calcFitness();
}

// This routine applies when -log2(rec) is integer 0,1,2,...
// rec = 1 => -log2 rec = 0 is OK here, each successive locus chosen from alternate parent
// assumes that random integer has random bits

void SetBabyGenotypeLogRec(Individual& Parent1, Individual& Parent2, Individual& baby)
{
    auto& g1 = Parent1.genotype;
    auto& g2 = Parent2.genotype;
    auto& gb = baby.genotype;
    auto& s1 = Parent1.stochast;
    auto& s2 = Parent2.stochast;
    auto& sb = baby.stochast;
    ulong rawint = rnd.rawint();
    ulong recShift = Individual::negLog2Rec; // -log 2 rec, w/rec = (1/2, 1/4, 1/8, ...), set in Popul
    ulong mask = (1 << recShift) - 1;        // e.g., recShift = 2 => mask = 00...0011, ie, low two bits
    ulong chrFlag = rawint & 1;              // determines initial parent w/prob = 1/2, ie, random bit
    auto rbits = rnd.bitSize() - recShift;   // remaining bits available
    
    for (int i = 0; i < Parent1.totalLoci; ++i){
        gb[i] = (chrFlag) ? g1[i] : g2[i];
        if (Individual::stoch) sb[i] = (chrFlag) ? s1[i] : s2[i];
        rawint >>= recShift;                            // move used bits out
        if ((rawint & mask) == mask) chrFlag ^= 1;      // flip flag if recombination
        if ((rbits -= recShift) == 0){                  // reload random bits if all used up
            rawint = rnd.rawint();                      // new random int
            rbits = rnd.bitSize();                      // reset remaining bits left to use
        }
    }
    baby.fitness = baby.calcFitness();
}

// No recombination, choose just one parent and copy genotype to baby. Maybe memcpy would be faster?? Note how to turn off warning for unused parameter

void SetBabyGenotypeNoRec(Individual& Parent, Individual& Unused __attribute__((unused)), Individual& baby)
{
    for (int i = 0; i < Parent.totalLoci; ++i){
        baby.genotype[i] = Parent.genotype[i];
        if (Individual::stoch) baby.stochast[i] = Parent.stochast[i];
    }
    baby.fitness = baby.calcFitness();
}

// Calculation of num and den take from openVclose.h in pagmo optimization code; assumes dentilde = den, ie, not studying role of variable plant w/regard to stability margin. Plant set, see manuscripts. Plant parameters do not vary, thus a is set to optimal value of a = sqrt(1 + gamma), and optimal value of J = sqrt(gamma).
// Forms for num and den in MMA file

double Individual::calcJ()
{
    double tmax = 20.0;
    double a = sqrt(1+gamma);
    if (abs(aSD) > 1e-6) a *= pow(2.0,rnd.normal(0,aSD));   // a = a*2^x, x ~ N(0,aSD)
    // p0 = 0 by assumption
    double p1 = genotype[0] * ((stoch) ? pow(2.0,rnd.normal(0,stochWt*stochast[0])) : 1.0);
    double p2 = genotype[1] * ((stoch) ? pow(2.0,rnd.normal(0,stochWt*stochast[1])) : 1.0);
    double q0 = genotype[2] * ((stoch) ? pow(2.0,rnd.normal(0,stochWt*stochast[2])) : 1.0);
    double q1 = genotype[3] * ((stoch) ? pow(2.0,rnd.normal(0,stochWt*stochast[3])) : 1.0);
    double q2 = genotype[4] * ((stoch) ? pow(2.0,rnd.normal(0,stochWt*stochast[4])) : 1.0);
    double r, k;
    if (loop == Loop::dclose){
        r = genotype[5] * ((stoch) ? pow(2.0,rnd.normal(0,stochWt*stochast[5])) : 1.0);
        k = genotype[6] * ((stoch) ? pow(2.0,rnd.normal(0,stochWt*stochast[6])) : 1.0);
    }
    std::vector<double> num;
    std::vector<double> den;
    switch (loop){
        case Loop::open:
            num = {q2,q1,q0};
            den = {p2, p1+a*p2, a*p1 + p2, p1};
            break;
        case Loop::close:
            num = {q2,q1,q0};
            den = {p2+q2, p1+a*p2+q1, a*p1+p2+q0, p1};
            break;
        case Loop::dclose:
            double rk = r*k;
            num = {rk*q2, rk*q1 + k*q2, rk*q0 + k*q1, k*q0};
            den = {rk*q2, p2 + rk*q1 + q2 + k*q2, p1 + a*p2 + rk*q0 + q1 + k*q1,
                a*p1 + p2 + q0 + k*q0, p1};
            break;
    }
    
    return performance(num, den, gamma, tmax, signalType::output);
}

double Individual::calcFitness()
{
    double optJ = sqrt(gamma);
    double Jdev = (calcJ()/optJ) - 1.0;
    return fitness = exp(-(Jdev*Jdev)/(2*fitVar));
}

