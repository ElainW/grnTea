import os
import global_values

class PrototypeParameters():
    def __init__(self):
        self.m_cov_par = {}  # each record in format (nclip, ndisc, nclip_disc)

    def get_par_by_cov(self, icov):
        i_cur_cov=int(icov)
        b_hit=False
        par_rcd=None
        for i in range(0,10):
            i_minus=i_cur_cov-i
            if i_minus in self.m_cov_par:
                par_rcd = self.m_cov_par[i_minus]
                b_hit=True
                break
            i_plus=i_cur_cov+i
            if i_plus in self.m_cov_par:
                par_rcd = self.m_cov_par[i_plus]
                b_hit = True
                break

        if b_hit==False:#no any hit
            i_tmp=int(icov)/10 * 10
            if i_tmp<10:
                return self.m_cov_par[5]
            elif i_tmp>200:
                return self.m_cov_par[200]
            else:
                return self.m_cov_par[100]
        return par_rcd

####
class Parameters(PrototypeParameters):
    def __init__(self):
        PrototypeParameters.__init__(self)
        #this is the recommended cutoff value based on coverage
        self.m_cov_par[5] = (1, 3, 0)
        self.m_cov_par[10]=(2,3, 0)
        self.m_cov_par[15] = (2, 3, 0)
        self.m_cov_par[20]=(2,4, 0)
        self.m_cov_par[25] = (2, 4, 1)
        self.m_cov_par[30]=(3,4, 1)
        self.m_cov_par[35] = (3, 5, 1)
        self.m_cov_par[40]=(3,6, 1)
        self.m_cov_par[45] = (4, 6, 1)
        self.m_cov_par[50] = (4, 7, 1)
        self.m_cov_par[60] = (5, 8, 2)
        self.m_cov_par[70] = (5, 9, 2)
        self.m_cov_par[80] = (6, 10, 2)
        self.m_cov_par[90] = (7, 11, 3)
        self.m_cov_par[100] = (8, 12, 3)
        self.m_cov_par[110] = (9, 12, 3)
        self.m_cov_par[120] = (9, 13, 3)
        self.m_cov_par[130] = (10, 14, 3)
        self.m_cov_par[140] = (11, 14, 4)
        self.m_cov_par[150] = (12, 15, 4)
        self.m_cov_par[175] = (13, 17, 5)
        self.m_cov_par[200] = (14, 20, 5)
        self.m_cov_par[250] = (20, 25, 7)
        self.m_cov_par[300] = (25, 30, 8)
########

