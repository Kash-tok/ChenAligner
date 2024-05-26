#include <iostream>
#include "ChenAligner.h"

constexpr char USE_SITUATION[] = "Fisrt Time";

int main()
{
    unsigned int alphabet_size = 5;  // A G C T !
    char erase_marker = '!';
    ChenAligner aligner(alphabet_size, erase_marker);

    if constexpr (USE_SITUATION == "First Time")
    {
        aligner.GetWatermarksFromFile("Your Watermark File.txt");
        aligner.DumpWatermarkStructureToCache("Where You Want to Store Your Structure.cache");
    }
    else
    {
        aligner.GetWatermarksFromCache("Where You Stored Your Structure.cache");
    }
    
    aligner.MatchReads(ReadType::FASTQ, "The Read File.fastq", "Where You Want to Store Results.csv", 0.5, 10);

    return 0;
}