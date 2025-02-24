import sys
import os
from circuit import *


class Code:

	stabilizer = [
                       ['X','Z','Z','X','I'],
                       ['I','X','Z','Z','X'],
                       ['X','I','X','Z','Z'],
                       ['Z','X','I','X','Z']]

	logical = {'X': ['X','X','X','X','X'],
		   'Z': ['Z','Z','Z','Z','Z']}


	stabilizer_CHP=[
			  '+XZZXI',
			  '+IXZZX',
			  '+XIXZZ',
			  '+ZXIXZ']


   

	destabilizer_CHP=[
			   '+ZIZII',
		   	   '+IIIIZ',
			   '+IIZII',
			   '+ZIZZI']
	
	stabilizer_logical_CHP = {'+Z': '+ZZZZZ',
				              '-Z': '-ZZZZZ'}

	destabilizer_logical_CHP = { '+Z':'+ZXZII',
				                 '-Z':'+ZXZII'}

	# destabs anticommute with all of the stabs except the corresponding stab for which anticommutes.
	# all stabs commute with themsevles and all destabs commute with themslves.

			



	# Option 3
	complete_look_up_table = { 
				 0: ['I','I','I','I','I'],
				 1: ['X','I','I','I','I'],
				 8: ['I','X','I','I','I'],	
				12: ['I','I','X','I','I'],
				 6: ['I','I','I','X','I'],
				 3: ['I','I','I','I','X'],
				11: ['Y','I','I','I','I'],
				13: ['I','Y','I','I','I'],
				14: ['I','I','Y','I','I'],
				15: ['I','I','I','Y','I'], 
			     7: ['I','I','I','I','Y'], 
	   			10: ['Z','I','I','I','I'],
				 5: ['I','Z','I','I','I'],	
				 2: ['I','I','Z','I','I'],
				 9: ['I','I','I','Z','I'], 
				 4: ['I','I','I','I','Z'] 	
				}  


	stabilizer_syndrome_dict = { 
				 tuple([0,0,0,0]): ['I','I','I','I','I'],
				 tuple([0,0,0,1]): ['X','I','I','I','I'],
				 tuple([1,0,0,0]): ['I','X','I','I','I'],	
				 tuple([1,1,0,0]): ['I','I','X','I','I'],
				 tuple([0,1,1,0]): ['I','I','I','X','I'],
				 tuple([0,0,1,1]): ['I','I','I','I','X'],
				 tuple([1,0,1,1]): ['Y','I','I','I','I'],
				 tuple([1,1,0,1]): ['I','Y','I','I','I'],
				 tuple([1,1,1,0]): ['I','I','Y','I','I'],
				 tuple([1,1,1,1]): ['I','I','I','Y','I'], 
			     tuple([0,1,1,1]): ['I','I','I','I','Y'], 
	   			 tuple([1,0,1,0]): ['Z','I','I','I','I'],
				 tuple([0,1,0,1]): ['I','Z','I','I','I'],	
				 tuple([0,0,1,0]): ['I','I','Z','I','I'],
				 tuple([1,0,0,1]): ['I','I','I','Z','I'], 
				 tuple([0,1,0,0]): ['I','I','I','I','Z'] 	
				}  





  
   
					
