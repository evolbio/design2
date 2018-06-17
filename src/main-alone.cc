#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <cctype>

#include <boost/filesystem.hpp>
#include <boost/asio/ip/host_name.hpp>

#include "fmt/format.h"

#include "param.h"
#include APPL_H

bool showProgress = false;
constexpr int maxLinesPerRun = 20;

/********************** Prototypes ****************************/

void 		    InitRuns(int& first, int& last, std::fstream& randFile, std::ifstream& paramFile,
                         std::ofstream& outFile, const std::string& exp);
void 		    UpdateRandFile(std::fstream& randFile, int first, int last);
std::string     WriteParmBuf(int first, std::fstream& randFile, std::ifstream& paramFile);

/**************************************************************/

int main(int argc, char *argv[])
{
	int i, first, last, arg;
    std::string exp;

    std::string usage =
        fmt::format("\n\tUSAGE:  {} -s experiment\n\n", argv[0])
        + "\t\t-s to show progress on stdout\n\n"
        + "\t\texperiment must begin with a letter\n\n";
    try {
        if (argc == 1) throw std::exception();
        // if '-' switch, assume -s and set arg = 2 to expect second arg
        if (argv[1][0] == '-'){
            showProgress = true;
            arg = 2;
        }
        // expect one arg, which is experiment letter
        else arg = 1;
        if (arg < argc && std::isalpha(argv[arg][0]))
            exp = argv[arg];
        else throw std::exception();
    }
    catch (const std::exception& e) {
        std::cerr << usage << std::endl;
        exit(1);
    }
    try {
        if (linesPerRun > maxLinesPerRun)
            ThrowError(__FILE__, __LINE__, "LinesPerRun greater than maxLinesPerRun.");
        MakeParam("input/", "design", exp.c_str(), 2);
        std::ifstream paramFile;
        std::fstream randFile;
        std::ofstream outFile;
        InitRuns(first, last, randFile, paramFile, outFile, exp);
        for (i = first; i <= last; ++i){
            std::istringstream parmBuf{WriteParmBuf(i, randFile, paramFile)};
            std::string result = Control(parmBuf);
            outFile << result;
            outFile.flush();
            UpdateRandFile(randFile, i+1, last);
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        std::ofstream("output/00_error.log") << e.what() << std::endl;
        exit(1);
    }
	return 0;
}

std::string WriteParmBuf(int first, std::fstream& randFile, std::ifstream& paramFile)
{
    randFile.seekg(0);
    std::string tmp;
    unsigned long seed;
    for (unsigned i = 0; i < 11; ++i) randFile >> tmp;
    randFile >> seed;
    std::string buf;
	for (int i = 0; i < linesPerRun; ++i){
        std::string lineBuf;
        std::getline(paramFile,lineBuf);
        buf += lineBuf;
	}
    std::string parmBuf = fmt::format("{} {} {} {}", first, first, seed, buf);
    return parmBuf;
}

void InitRuns(int& first, int& last, std::fstream& randFile, std::ifstream& paramFile,
              std::ofstream& outFile, const std::string& exp)
{
    auto hostname_fqdn = boost::asio::ip::host_name();
    auto hostname_short = hostname_fqdn.substr(0,hostname_fqdn.find_first_of('.'));
    // might want to split short name on '-' for some clusters that use nXX-host names

    std::string filename = fmt::format("{}.{}", "input/random", hostname_short);
    randFile.open(filename);
    if (!randFile)
        ThrowError(__FILE__, __LINE__, "Could not open " + filename);
    std::string tmp;
    randFile >> tmp >> tmp >> tmp >> first;
    randFile >> tmp >> tmp >> tmp >> last;
	if (randFile.bad())
        ThrowError(__FILE__, __LINE__, "Error reading iterates for main control loop from " + filename);
	if (first <= 0)
        ThrowError(__FILE__, __LINE__, "First run should be 1 or greater, do not use 0 in " + filename);
    filename = fmt::format("input/Exp{}.parm.{}", exp, hostname_fqdn);
    std::string filename2 = fmt::format("input/Exp{}.parm.{}", exp, hostname_short);
    boost::filesystem::rename(filename, filename2);
	paramFile.open(filename2);
    for (int i = 0; i < linesPerRun * first; ++i)
        std::getline(paramFile,tmp);
    filename = fmt::format("output/data.Exp{}.{}", exp, hostname_short);
    if (boost::filesystem::exists(filename)) {
        filename2 = fmt::format("{}.bak", filename);
        boost::filesystem::rename(filename, filename2);
    }
	outFile.open(filename);
}

void UpdateRandFile(std::fstream& randFile, int first, int last)
{
    randFile.seekp(0);
    randFile << fmt::format("Current run = {}\n", first);
    randFile << fmt::format("Last run    = {}\n", last);
    randFile << fmt::format("Rand seed   = {}\n{:20}", rnd.rawint(), "");
    randFile.flush();
}


