#! /bin/sh
#
# Copyright (C) 2020 Chedly Kastally <ckastall@gmail.com>
# Distributed under terms of the MIT license.
#
# Handle the output of paralog_window_filtering.R to produce a bed file.

window_half_size=125
window_threshold=0.9
# For verbose output, set value to 1
verbose_mode=0

awk -v verbose="${verbose_mode}" \
    -v whalfsize="${window_half_size}" \
    -v thresh="${window_threshold}" '
BEGIN{
    FS="\t"
    OFS="\t"
    OFMT="%.0f" # this will round the position value
    score_current=""
    pos_current=""
    contig_current=""
    filter_current=""
}
(NR==1){
    if($3  != "contig" ||
       $4  != "pos"    ||
       $17 != "filter" ||
       $18 != "window_score_proportion_good_snp" ||
       $19 != "window_filter_TrueIsRetained"){
        print "ERROR: Format error, header is not what it should be. Aborting." ;
        exit 1 ;
    } else {
        if (verbose) {
            print "# File format is correct; moving on." ; 
        }
    }}
(NR==2) {
    contig_current=$3
    pos_current=$4
    status_current=$17
    score_current=$18
    filter_current=$19
    current_state_window = $18 > thresh ? "TRUE" : "FALSE"
    win_start=($4 - 1 - whalfsize)
    win_end=($4 + whalfsize)
}
(NR>2) {

    if (verbose) {

        print "# examining position:", $4 ;
        print "# contig_current:", contig_current ;
        print "# pos_current", pos_current ;
        print "# current window start:", win_start ;
        print "# current window end:", win_end ;

    }

    if ($17 == "unknown") { next ;}

    if ($4 > pos_current + (whalfsize*2) ||
        $3 != contig_current) {
        
        # Moving away from the current window: print it and move on

        if (current_state_window == "FALSE") {
            if (win_start < 1) {win_start = 1}
            print contig_current, win_start, win_end
        }

        if (verbose) {
            print "# initialzing new window, at contig and pos", $3, $4 ;
        }

        win_start=($4 - 1 - whalfsize)
        win_end=($4 + whalfsize)

    } else {

        # We are still in a window being evaluated, but as we move, the score
        # change depending on surrounding SNPs

        state_curr_position_in_window = $18 > thresh ? "TRUE" : "FALSE"

        if (state_curr_position_in_window != current_state_window) {

            # Changing the status of the current window.

            if (verbose) { print "# Transition" ; }
            if (state_curr_position_in_window == "FALSE") {

                if (verbose) { print "# Going into an excluded area" ; }

                # If going from TRUE to FALSE, we go from included into an excluded area.
                # The included area does not need printing, we can start a new
                # window adjust the start position
                #
                # Here we define the new starting point. This should be negative
                # the formula is basically a difference on the y-axis divided by the slope
                win_start_adjust = (($18 - thresh) / ((score_current - $18) / (pos_current - $4)))

                # This should be negative
                if (win_start_adjust < 0) {
                    print "there is a problem here:", pos_current, win_start_adjust;
                    exit 1;
                }

                # Cap the adjustment to -whalfsize, ie, we do not exclude areas
                # of the genome more than whalfsize from an SNP
                if (win_start_adjust > whalfsize) {win_start_adjust = whalfsize}

                # -1 because bed files are 0-padded
                win_start = ($4 - win_start_adjust) - 1

                # Edge case: the correction is such that the previous position
                # (good), is now in the excluded area (as if it were now bad).
                #
                # This happens notably if you go from a score of > thresh to < thresh
                # within a very short distance, eg 1 bp;
                #
                # for those, I fix the limit to right after the last good position (pos_current - 1).
                # Another option would be to take the middle point..

                if (win_start < pos_current) {

                    # print "edgecase!" ;

                    win_start = pos_current
                }

                win_end=($4 + whalfsize)

                if (verbose) {
                    ## Debug print; to check the coordinate calculation of the new starting position of the window
                    print "# Adjusting starting position of window:";
                    print "# thresh:", thresh ;
                    print "# Ya, score current:", score_current ;
                    print "# Yb, win score: ", $18 ;
                    print "# xa, position current:", pos_current ;
                    print "# xb, position line:", $4 ;
                    print "# numera", (score_current - $18) ;
                    print "# denomin", (pos_current - $4) ;
                    print "# start adjustment: ", win_start_adjust ;
                    print "# start after adjustment: ", win_start ;
                }

            } else {

                # From FALSE to TRUE,
                # end a window in transition region, adjust the end transition
                if (verbose) { print "# Leaving an excluded area" ; }

                # Define the new ending point.
                # do note that there is a small difference in the formula compared to the symetrical case: it is normal, we do not measure exactly the same thing
                win_end_adjust = ((thresh - score_current) / (($18 - score_current) / ($4 - pos_current)))

                if (win_end_adjust < 0) {
                    print "there is a problem here:", pos_current, win_end_adjust;
                    exit 1;
                }

                # Cap the adjustment to +whalfsize
                if (win_end_adjust > whalfsize) {win_end_adjust = whalfsize}

                win_end = (pos_current + win_end_adjust)

                # Edge case: the correction is such that the line position
                # (good), is now in the excluded area (as if it were now bad).
                #
                # This happens notably if you go from a score of < thresh to > thresh
                # within a very short distance, eg 1 bp;
                #
                # for those, I fix the limit to right before the next good position ($4 - 1).
                # Another option would be to take the middle point.

                if (win_end > ($4 - 1)) {

                    # print "edgecase!" ;

                    win_end = $4 - 1
                }

                if (verbose) {
                    ## Debug print; to check the coordinate calculation of the new ending position of the window
                    print "# Adjusting ending position of window:" ;
                    print "# thresh:", thresh ;
                    print "# Ya, score current:", score_current ;
                    print "# Yb, win score: ", $18 ;
                    print "# xa, position current:", pos_current ;
                    print "# xb, position line:", $4 ;
                    print "# numerator", (score_current - $18) ;
                    print "# denominator", (pos_current - $4) ;
                    print "# end adjustment: ", win_end_adjust ;
                    print "# end after adjustment: ", win_end ;
                }

                if (win_start < 1) {win_start = 1}
                print contig_current, win_start, win_end

                win_start=win_end

            }
        } else {

            win_end=($4 + whalfsize)
        }
    }

    if ($17 == "excluded" &&
        $18 > thresh) {

        if (verbose) {
            print "# Bad SNP in good window detected; contig, position:", $3, $4 ;
        }

        print $3, $4 - 1, $4 ;
    }

    contig_current=$3
    pos_current=$4
    score_current=$18
    filter_current=$19
    current_state_window = $18 > thresh ? "TRUE" : "FALSE"

}
END{
    if (current_state_window == "FALSE") {
        if (win_start < 1) {win_start = 1}
        print contig_current, win_start, win_end
    }}'
