# This module designed for testing functions in the command line


import utility
rp_pairs = [('Y',-1),('F',0), ('A', 3)]
rexp = '.F.A.{2,9}Y'
rp_ids = utility.rp_iteration_id(rp_pairs, rexp)

print "rp_ids=", rp_ids


