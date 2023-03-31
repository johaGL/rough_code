

lili = scipy.stats.mannwhitneyu(np.array([11, 15, 17, 18]), np.array([20, 24, 30, 35]),
                            use_continuity=False,
                            alternative="two-sided",
                            )
                            
# MannwhitneyuResult(statistic=0.0, pvalue=0.02857142857142857)                            

#  now a subset of each group 
lilisu = scipy.stats.mannwhitneyu(np.array([17, 18]), np.array([ 30, 35]),
                            use_continuity=False,
                            alternative="two-sided",
                            )
                            
# MannwhitneyuResult(statistic=0.0, pvalue=0.3333333333333333)



# ok for at least 3 samples:
lo =  scipy.stats.mannwhitneyu(np.array([14, 16, 18]), np.array([ 15, 17, 18]),
                            use_continuity=False,
                            alternative="two-sided",
                            )
# MannwhitneyuResult(statistic=3.5, pvalue=0.6579050194284821)



