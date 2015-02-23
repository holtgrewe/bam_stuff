// ==========================================================================
//                            BAM Duplicate Masker
// ==========================================================================
// Copyright (c) 2015, Manuel Holtgrewe, Charite Universitaetsmedizin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@charite.de>
// ==========================================================================

#include "bam_markdup_options.h"

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>

#include "version.h"

// ----------------------------------------------------------------------------
// Class BamMarkDupOptions
// ----------------------------------------------------------------------------

void BamMarkDupOptions::print(std::ostream & out) const
{
    out << "__OPTIONS________________________________________________________\n"
        << "\n"
        << "VERBOSITY       \t" << verbosity << "\n"
        << "\n"
        << "INPUT BAM   \t" << inPath << "\n"
        << "OUTPUT BAM  \t" << outPath << "\n"
        << "WINDOW SIZE \t" << windowLength << "\n";
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

BamMarkDupOptions parseCommandLine(int argc, char ** argv)
{
    BamMarkDupOptions result;

	for (int i = 0; i < argc; ++i)
		result.commandLine.push_back(argv[i]);

    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bam_markdup");

    // Set short description, version, and date.
    setShortDescription(parser, "BAM Duplicate Marker");
    setVersion(parser, PROGRAM_VERSION);
    setDate(parser, PROGRAM_DATE);

    // Define usage line and long description.
    addUsageLine(parser, "--in IN.bam --out OUT.bam");
    addDescription(parser, "Mark duplicates in IN.bam to OUT.bam");

    addOption(parser, seqan::ArgParseOption("q",  "quiet",        "Quiet output"));
    addOption(parser, seqan::ArgParseOption("v",  "verbose",      "Verbose output"));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose output"));

    // Define Options -- Section Input / Output Optiosn
    addSection(parser, "Input / Output Options");

    addOption(parser, seqan::ArgParseOption("", "in", "Input BAM alignment.",
                                            seqan::ArgParseArgument::INPUT_FILE, "BAM"));
    setRequired(parser, "in", true);
    setValidValues(parser, "in", seqan::BamFileIn::getFileExtensions());

    addOption(parser, seqan::ArgParseOption("", "out", "Output BAM file.",
                                            seqan::ArgParseArgument::INPUT_FILE, "BAM"));
    setRequired(parser, "out", true);
    setValidValues(parser, "out", seqan::BamFileOut::getFileExtensions());

    // Define Options -- Algorithm Parameters
    addSection(parser, "Algorithm Parameters");

    addOption(parser, seqan::ArgParseOption("", "window-length", "Window length.",
                                            seqan::ArgParseArgument::INTEGER, "LEN"));
    setDefaultValue(parser, "window-length", 4000);

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        throw InvalidCommandLineArgumentsException();

    // Extract option values.
    result.verbosity = isSet(parser, "quiet") ? 0 : result.verbosity;
    result.verbosity = isSet(parser, "verbose") ? 2 : result.verbosity;
    result.verbosity = isSet(parser, "very-verbose") ? 3 : result.verbosity;

    getOptionValue(result.inPath, parser, "in");
    getOptionValue(result.outPath, parser, "out");

    getOptionValue(result.windowLength, parser, "window-length");

    return result;
}
