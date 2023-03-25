import numpy as np

def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def hypergeom_enrichment_score(MIA_qvalue_list):
    MIA_qvalue_list_log = []
    for i in range(len(MIA_qvalue_list)):
        enr = -np.log10(MIA_qvalue_list[i])
        dep = -np.log10(1 - MIA_qvalue_list[i])
        # print(enr,dep)
        if enr < dep:
            MIA_qvalue_list_log.append(-dep)
        else:
            MIA_qvalue_list_log.append(enr)

    # Replace infinite with max/min if exists
    enr_max = np.max([i for i in MIA_qvalue_list_log if np.isfinite(i)])
    enr_min = np.min([i for i in MIA_qvalue_list_log if np.isfinite(i)])
    MIA_qvalue_list_log = [i if np.isfinite(i) else enr_max if i > 0 else enr_min for i in MIA_qvalue_list_log]
    return MIA_qvalue_list_log