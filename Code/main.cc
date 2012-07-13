// =====================================================================================
// 
//       Filename:  main.cc
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  04/06/2011 19:49:32
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Jianxing Feng (feeldead), feeldead@gmail.com
//        Company:  Tongji Univ.
// 
// =====================================================================================

#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include "GeneInfo.hpp"
#include "GFOLD.hpp"
#include "Utility.hpp"

#define VERSION "V1.0.4"
#define DATE "Sat Apr 28 16:01:50 CST 2012"


using namespace std;

void
Help()
{
    cout << endl;
    cout << "      =============================================================================== " << endl;
    cout << "          gfold    :   Generalized fold change for ranking differentially expressed   " << endl;
    cout << "                       genes from RNA-seq data." << endl;
    cout <<                                                     endl;
    cout << "          Author   :   Jianxing Feng (jianxing.tongji@gmail.com)" << endl;
    cout << "            Date   :   " << DATE << endl; 
    cout << "          Version  :   " << VERSION << endl;
    cout << "      =============================================================================== " << endl;
    cout <<                                                   endl;
    cout << "      USAGE:   Please use command 'man doc/gfold.man' to find documentation." << endl;
    cout <<                                endl;
}

void
die_arg_missing(bool b_die, string argument)
{
    if (b_die)
    {
        cerr << "ERROR: argument " << argument << " is required" << endl;
        Help();
        exit(0);
    }   
}

void
die_arg_wrong(bool b_die, string argument)
{
    if (b_die)
    {
        cerr << "ERROR: Parameter for argument " << argument << " is invalid. Please refer to the manual" << endl;
        Help();
        exit(0);
    }   
}


int 
main(int argc, char* argv[])
{
    int startArg = 1;

    string job = "";

    // For job 'preprocess'
    string ref_seq_file = "";
    string gene_annotation = "";
    string gene_annotation_format = "GPF";
    string short_reads_file = ""; 
    string short_reads_format = "SAM"; 

    // For job 'estimate'
    int strand_specific_code = 0;
    string output_file = "";
    string sample_suffix = "";

    int verbos_level = 2;
    int burn_in_count = 1000;
    int sampled_count = 1000;
    double significant_cutoff = 0.01;        
    int random_sampled_pairs = 20;
    bool b_accurate = true;

    // For job 'diff'
    vector<string> first_group_samples;
    vector<string> second_group_samples;
    string normal_method = "DESeq";

    // Begin parsing the jobs and parameters
    for (int i = startArg; i < argc; i++)
    {
        if (strcmp(argv[i], "-h") == 0)
        {
            Help();
            return 0;
        }

        else if (strcmp(argv[i], "count") == 0)
            job = "count";
        else if (strcmp(argv[i], "diff") == 0)
            job = "diff";
        else if (strcmp(argv[i], "-ann") == 0)
            gene_annotation = argv[++i];
        else if (strcmp(argv[i], "-annf") == 0)
            gene_annotation_format = argv[++i];
        else if (strcmp(argv[i], "-tag") == 0)
            short_reads_file = argv[++i];
        else if (strcmp(argv[i], "-tagf") == 0)
            short_reads_format = argv[++i];
        else if (strcmp(argv[i], "-o") == 0)
            output_file = argv[++i];
        else if (strcmp(argv[i], "-suf") == 0)
            sample_suffix = argv[++i];
        else if (strcmp(argv[i], "-s") == 0)
            if (argv[++i][0] == 'T') strand_specific_code = 1;
            else strand_specific_code = 2;
        else if (strcmp(argv[i], "-acc") == 0)
	    b_accurate = (argv[++i][0] == 'T');
        else if (strcmp(argv[i], "-si") == 0)
            sampled_count = atoi(argv[++i]);
        else if (strcmp(argv[i], "-bi") == 0)
            burn_in_count = atoi(argv[++i]);
        else if (strcmp(argv[i], "-v") == 0)
            verbos_level = atoi(argv[++i]);
        else if (strcmp(argv[i], "-s1") == 0)
        {
            vector<string> temp;
            split(argv[++i], ',', temp);
            first_group_samples.insert(first_group_samples.end(), temp.begin(), temp.end());
        }
        else if (strcmp(argv[i], "-s2") == 0)
        {
            vector<string> temp;
            split(argv[++i], ',', temp);
            second_group_samples.insert(second_group_samples.end(), temp.begin(), temp.end());
        }
        else if (strcmp(argv[i], "-norm") == 0)
            normal_method = argv[++i];
        else if (strcmp(argv[i], "-sc") == 0)
            significant_cutoff = atof(argv[++i]);
        else if (strcmp(argv[i], "-r") == 0)
            random_sampled_pairs = atoi(argv[++i]);
        else 
        {
            Help();
            cerr << "Wrong parameter " << argv[i] << endl;
            exit(0);
        }
    }

    die_arg_missing(job == "", "Job"); 
    die_arg_missing(output_file == "", "-o"); 

    if (job == "count")
    {
        die_arg_missing(gene_annotation == "", "-ann"); 
        die_arg_missing(short_reads_file == "", "-tag"); 
        die_arg_wrong(gene_annotation_format != "GPF" && gene_annotation_format != "BED", "-annf");
        die_arg_wrong(short_reads_format != "SAM" && short_reads_format != "BED", "-tagf");

        GeneInfo gene_info(output_file, 1, strand_specific_code == 1, false, verbos_level);
        if (gene_annotation_format == "GPF")
            scanGPF(gene_annotation, gene_info, verbos_level);
        else
            scanAnnotBED(gene_annotation, gene_info, verbos_level);

        if (short_reads_format == "SAM")
            scanSAM(short_reads_file, gene_info, verbos_level);
        else
            scanReadsBED(short_reads_file, gene_info, verbos_level);

        gene_info.PrintGeneReadCount();
    }
    else if (job == "diff")
    {
        die_arg_missing(first_group_samples.size() == 0, "-s1"); 
        die_arg_missing(second_group_samples.size() == 0, "-s2"); 

        if (normal_method != "Count" && normal_method != "TMM" && normal_method != "DESeq")
        {
            cerr << "ERROR: Unknown normalization method : " << normal_method << endl;
            exit(1);
        }

        if (first_group_samples == second_group_samples)
        {
            cerr << "ERROR: Two input groups of samples are exactly the same!" << endl;
            exit(1);
        }

        GFOLD gfold(VERSION, verbos_level, normal_method, burn_in_count, sampled_count, significant_cutoff, random_sampled_pairs, b_accurate);
        gfold.CalculateAll(first_group_samples, second_group_samples, sample_suffix, output_file);
    }
    cerr << "Job " << job << " is DONE!" << endl;

    return 0;
}
