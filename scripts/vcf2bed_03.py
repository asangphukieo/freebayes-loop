import sys
import gzip

def find_deletions(min_del, vcf_input, split_thres):
    """
    Finds deletions in a VCF file that exceed a minimum deletion length.

    Args:
        min_del (int): Minimum deletion length (bp) to print out.
        vcf_input (str): Path to the VCF input file.
        split_thres (int): Threshold in bp to split large deletions into chunks.

    Returns:
        list: List of tuples containing chromosome, start position, and end position of deletions.

    Process:
        1. It iterates over each line in the VCF file using a for loop.
        2. If the line does not start with a "#" (indicating a header line), it proceeds with further processing.
        3. The line is split into columns using the tab ('\\t') delimiter, and the resulting columns are stored in the col variable.
        4. The length of the reference allele (col[3]) is calculated and stored in the del_len variable.
        5. The script checks if the deletion length (del_len) is greater than the specified minimum deletion length (min_del) and if the seventh column (INFO field) contains an item with "TYPE=del" indicating a deletion.
        6. If both conditions are met, the script extracts the chromosome (col[0]), start position (col[1]), and calculates the end position (start position + deletion length). It then prints the chromosome, start position, and end position.
        7. The position is collected in deletions list , if it not been collected before.
        8. check position of the next line to have no longer than min_del , if it is longer than min_del, it may be large deletion

    Version:
        03 : - Add option to separate deletion chunk if it larger than ... bp as user input option;
             - edit end position to end_pos - 1 (bug)
             - edit the starting point of each chunk -1, because freebayes does not calculate variant for first position
        02 : Remove option to check TYPE= deletion to include every large variant, to avoid variant calling skipping
             Add option to separate deletion chunk if it larger than 100 bp
    """
    deletions = []

    #vcf file need to be sorted !
    with gzip.open(vcf_input,"rt") as file:
        line = file.read().strip().split('\n')
        n=0
        while n < len(line):
            if n == len(line)-1 : #if it is last line
                if not line[n].startswith("#") and next_line != "" :
                    col = line[n].split('\t')
                    del_len = len(col[3])
                    if del_len > min_del :
                        chrom = col[0]
                        start_pos = int(col[1])
                        end_pos = start_pos + del_len

                        if del_len > split_thres : #split chunks
                            ori_start = start_pos
                            for no in range(1,int(del_len/split_thres)+1):
                                if no == 1 : #first chunk
                                    start_pos = start_pos
                                    end_pos = ori_start + split_thres  #-1
                                elif no != int(del_len/split_thres):
                                    start_pos = end_pos #+ 1
                                    end_pos = start_pos + split_thres #-1
                                elif no == int(del_len/split_thres):
                                    if (int(del_len % split_thres) ==0) :
                                        start_pos = end_pos #+ 1
                                        end_pos = start_pos + split_thres
                                    else:
                                        start_pos = end_pos #+ 1
                                        end_pos = start_pos + split_thres #-1

                                if (chrom, start_pos, end_pos) not in deletions: #check duplicate
                                    deletions.append((chrom, start_pos, end_pos))

                                if no == int(del_len/split_thres) : #in case of last chunk
                                    if (int(del_len % split_thres) >=1) : #in case of last chunk and have remainder
                                        start_pos = end_pos #+ 1
                                        end_pos = ori_start + ( (no * split_thres)+ int(del_len % split_thres))

                                    if (chrom, start_pos, end_pos) not in deletions: #check duplicate
                                        deletions.append((chrom, start_pos, end_pos))
                            else:
                                if (chrom, start_pos, end_pos) not in deletions: #check duplicate
                                        deletions.append((chrom, start_pos, end_pos))
            else:
                next_line = line[n+1]
                if not line[n].startswith("#") and next_line != "" :
                    col = line[n].split('\t')
                    del_len = len(col[3])
                    if del_len > min_del :
                            if (int(next_line.split()[1]) >= int(line[n].split()[1]) + min_del) :  #if next line does not in range of min_del
                                chrom = col[0]
                                start_pos = int(col[1])
                                end_pos = start_pos + del_len

                                if del_len > split_thres : #split chunks
                                    ori_start = start_pos
                                    for no in range(1,int(del_len/split_thres)+1):
                                        if no == 1 : #first chunk
                                            start_pos = start_pos
                                            end_pos = ori_start + split_thres  #-1
                                        elif no != int(del_len/split_thres):
                                            start_pos = end_pos #+ 1
                                            end_pos = start_pos + split_thres #-1
                                        elif no == int(del_len/split_thres):
                                            if (int(del_len % split_thres) ==0) :
                                                start_pos = end_pos #+ 1
                                                end_pos = start_pos + split_thres
                                            else:
                                                start_pos = end_pos #+ 1
                                                end_pos = start_pos + split_thres #-1

                                        if (chrom, start_pos, end_pos) not in deletions: #check duplicate
                                            deletions.append((chrom, start_pos, end_pos))

                                        if no == int(del_len/split_thres) : #in case of last chunk
                                            if (int(del_len % split_thres) >=1) : #in case of last chunk and have remainder
                                                start_pos = end_pos #+ 1
                                                end_pos = ori_start + ( (no * split_thres)+ int(del_len % split_thres))

                                            if (chrom, start_pos, end_pos) not in deletions: #check duplicate
                                                deletions.append((chrom, start_pos, end_pos))

                                else:
                                    if (chrom, start_pos, end_pos) not in deletions: #check duplicate
                                            deletions.append((chrom, start_pos, end_pos))

            n=n+1

    return deletions

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python vcf2bed_03.py <min_del> <vcf_input> <bp_per_chunk>")
        sys.exit(1)

    min_del = int(sys.argv[1])
    vcf_input = sys.argv[2]
    split_thres = int(sys.argv[3])

    results = find_deletions(min_del, vcf_input, split_thres)
    if len(results) != 0 :
        for chrom, start_pos, end_pos in results:
                print(chrom, start_pos, end_pos)
    else:
        print("", end='')
