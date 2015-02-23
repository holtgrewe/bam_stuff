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

#include "bam_markdup_app.h"

#include <iostream>
#include <vector>
#include <string>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/simple_intervals_io.h>

#include "version.h"
#include "bam_markdup_options.h"

namespace {  // anonymous namespace

// Original comparison function for strings by Heng Li from samtools.

inline int strnum_cmp(const char *a, const char *b)
{
    char *pa, *pb;
    pa = (char*)a; pb = (char*)b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            long ai, bi;
            ai = strtol(pa, &pa, 10);
            bi = strtol(pb, &pb, 10);
            if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
        } else {
            if (*pa != *pb) break;
            ++pa; ++pb;
        }
    }
    if (*pa == *pb)
        return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;
    return *pa<*pb? -1 : *pa>*pb? 1 : 0;
}

bool isLeftOf(seqan::BamAlignmentRecord const & lhs, seqan::BamAlignmentRecord const & rhs)
{
    return ((unsigned)lhs.rID < (unsigned)rhs.rID ||
            (lhs.rID == rhs.rID && lhs.beginPos < rhs.beginPos) ||
            (lhs.rID == rhs.rID && lhs.beginPos == rhs.beginPos && strnum_cmp(toCString(lhs.qName), toCString(rhs.qName)) < 0));
}

}  // anonymous namespace

// ---------------------------------------------------------------------------
// Class BamMarkDupAppImpl
// ---------------------------------------------------------------------------

class BamMarkDupAppImpl
{
public:
    BamMarkDupAppImpl(BamMarkDupOptions const & options) :
        options(options), bamFileIn(bamFileOut)
    {}

    void run();

private:

    // Process regions one-by-one.
    void processRecords();

    // Open files.
    void openFiles();

    // Program configuration.
    BamMarkDupOptions options;

    // Objects used for I/O.
    seqan::BamFileIn bamFileIn;
    seqan::BamFileOut bamFileOut;

    // BAM header is read into this variable.
    seqan::BamHeader bamHeader;
};

void BamMarkDupAppImpl::run()
{
    double startTime = seqan::sysTime();

    if (options.verbosity >= 1)
    {
        std::cerr << "BAM Duplicate Masker\n"
                  << "====================\n\n";
        options.print(std::cerr);
    }

    // Open Files
    openFiles();

    // Process BAM Records

    if (options.verbosity >= 1)
        std::cerr << "\n"
                  << "__PROCESSING RECORDS_____________________________________________\n"
                  << "\n";

    processRecords();

    if (options.verbosity >= 1)
        std::cerr << "\nTotal time: " << (seqan::sysTime() - startTime) << " s\n"
                  << "All done. Have a nice day!\n";
}

void BamMarkDupAppImpl::openFiles()
{
    if (options.verbosity >= 1)
    {
        std::cerr << "\n"
                  << "__OPENING FILES__________________________________________________\n"
                  << "\n";
    }

    std::cerr << "Opening " << options.inPath << " ...\n";
    if (!open(bamFileIn, options.inPath.c_str()))
        throw new std::runtime_error("Could not open input file.");
    std::cerr << " OK\n";

    readHeader(bamHeader, bamFileIn);
    seqan::BamHeaderRecord pgRecord;
    pgRecord.type = seqan::BAM_HEADER_PROGRAM;
    appendValue(pgRecord.tags, seqan::BamHeaderRecord::TTag("ID", "bam_markdup"));
    appendValue(pgRecord.tags, seqan::BamHeaderRecord::TTag("PN", "bam_markdup"));
    std::string commandLine;
    for (auto token : options.commandLine)
    {
        if (!commandLine.empty())
            commandLine.append(" ");
        commandLine.append(token);
    }
    appendValue(pgRecord.tags, seqan::BamHeaderRecord::TTag("VN", PROGRAM_VERSION));
    appendValue(pgRecord.tags, seqan::BamHeaderRecord::TTag("CL", commandLine.c_str()));
    appendValue(pgRecord.tags, seqan::BamHeaderRecord::TTag("DS", "masking of duplicate alignments"));

    std::cerr << "Opening " << options.outPath << " ...\n";
    if (!open(bamFileOut, options.outPath.c_str()))
        throw new std::runtime_error("Could not open output file.");
    std::cerr << " OK\n";
    writeHeader(bamFileOut, bamHeader);
}

void BamMarkDupAppImpl::processRecords()
{
    double startTime = seqan::sysTime();



    if (options.verbosity >= 1)
        std::cerr << "Masking time: " << (seqan::sysTime() - startTime) << " s\n";
}

// ---------------------------------------------------------------------------
// Class BamMarkDupApp
// ---------------------------------------------------------------------------

BamMarkDupApp::BamMarkDupApp(BamMarkDupOptions const & options) :
        impl(new BamMarkDupAppImpl(options))
{}

BamMarkDupApp::~BamMarkDupApp()  // for pimpl
{}

void BamMarkDupApp::run()
{
    impl->run();
}
