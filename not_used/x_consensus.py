##11/22/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

# YW 2020/07/05 note: define python class ClusterChecker
# YW 2020/07/14 note: added + 1 to the RHS of _is_disc_cluster function: peak_win_cnt = i_max_right - i_max_left

class ClusterChecker():
    ####
    def _is_disc_cluster(self, l_pos, i_dist, ratio):
        l_pos.sort()  ###sort the candidate sites
        n_all_cnt = len(l_pos)
        if n_all_cnt == 0:
            return False, -1, -1
        elif n_all_cnt == 1: # YW 2020/07/14: to deal with the case when there is one element in l_pos
            return True, 0, 0
        i_max_left, i_max_right = self.find_max_cover(l_pos, i_dist)
        peak_win_cnt = i_max_right - i_max_left + 1 # +1 added by YW 2020/07/14, otherwise, if there's only one read, it is filtered out
        l_boundary = l_pos[i_max_left]
        r_boundary = l_pos[i_max_right]

        n_hit = peak_win_cnt
        if (float(n_hit) / float(n_all_cnt)) < ratio:
            return False, -1, -1
        return True, l_boundary, r_boundary
    ####
    #most of the left disc should be non_rc, while most of the right disc should be rc
    def _is_disc_orientation_consistency(self, n_l_rc, n_l_non_rc, n_r_rc, n_r_non_rc, ratio):
        b_l_consist=True
        if ((n_l_rc+n_l_non_rc) > 0) and (float(n_l_non_rc)/float(n_l_rc+n_l_non_rc) < ratio):
            b_l_consist = False
        b_r_consist = True
        if ((n_r_rc + n_r_non_rc) > 0) and (float(n_r_rc) / float(n_r_rc + n_r_non_rc) < ratio):
            b_r_consist = False
        return b_l_consist and b_r_consist

    ####
    # YW 2020/07/14 clarified this comment: return value is the longest interval of l_pos with l_pos[i_end]-l_pos[i_start]<=win_size (if there are >1 such intervals with the same longest length, the first will show up)
    # YW 2020/07/14 added more comments for the variables
    def find_max_cover(self, l_pos, win_size):
        i_right = 0
        i_total = len(l_pos)
        i_cnt = 0 # counts i_right - i_left with l_pos[i_right]-l_pos[i_left]<=win_size
        i_max = 0 # max of i_cnt, this is to help find the longest sub-list with l_pos[i_end]-l_pos[i_start]<=win_size
        i_start = 0
        i_end = 0
        for i_left in range(i_total):
            cur_val = l_pos[i_left]
            if i_right < i_left:
                i_right = i_left
            for i_temp in range(i_right, i_total):
                right_val = l_pos[i_temp]
                i_right = i_temp
                if (right_val - cur_val) > win_size:
                    if i_max < i_cnt:
                        i_max = i_cnt
                        i_start = i_left
                        i_end = i_right - 1
                    break
                i_cnt += 1
            if i_max < i_cnt:  # for the last one
                i_max = i_cnt
                i_start = i_left
                i_end = i_right

            i_cnt -= 1
            if i_cnt < 0:
                i_cnt = 0
            if i_right == (i_total - 1):
                break
        return i_start, i_end

    ##find the first and second peak clip position
    # ins_chrm is the ori reported insertion chrm on the reference
    # ins_pos is the ori reported insertion pos on the reference
    def find_first_second_peak(self, m_clip_pos, ins_chrm, ins_pos, BIN_SIZE):
        ##first check whether m_clip_pos[ins_chrm][ins_pos] exist
        peak_pos = -1
        peak_acm = 0
        scnd_peak_pos = -1
        scnd_peak_acm = 0
        mpos = m_clip_pos[ins_chrm][ins_pos]  # is a dictionary {ref_pos:[(cns_pos)]}
        istart = -1 * (BIN_SIZE / 2)
        iend = BIN_SIZE / 2
        for pos in mpos:  # cns_pos may be -1
            if pos <= 0:
                continue
            cur_acm = 0
            for offset in range(istart, iend):
                tmp_pos = pos + offset
                if tmp_pos in mpos:
                    cur_acm += len(mpos[tmp_pos])  ###number of reads clipped at the ref position
            if cur_acm > peak_acm:
                scnd_peak_pos = peak_pos
                scnd_peak_acm = peak_acm
                peak_acm = cur_acm
                peak_pos = pos
            elif cur_acm > scnd_peak_acm:
                scnd_peak_acm = cur_acm
                scnd_peak_pos = pos

        l_peak_nbr = []
        l_scnd_peak_nbr = []
        for offset in range(istart, iend):
            tmp_pos = peak_pos + offset
            if tmp_pos in mpos:
                l_peak_nbr.append(tmp_pos)

            tmp_scnd_pos = scnd_peak_pos + offset
            if tmp_scnd_pos in mpos:
                l_scnd_peak_nbr.append(tmp_scnd_pos)

        return peak_pos, l_peak_nbr, peak_acm, scnd_peak_pos, l_scnd_peak_nbr, scnd_peak_acm
