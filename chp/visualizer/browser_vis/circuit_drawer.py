#!/usr/bin/python

# import cgi
# import cgitb
#cgitb.enable()
from . import circuit
from . import gate_drawer
from .htmler import *

size = 18
boxx = size
boxy = size
separationy = size*1.5
separationx = size*2//3
linewidth = 1
qid_to_y = {}
time_to_x = {}

linkmax = 0
links = {}
qsize = 0
tsize = 0
MAXDATA = 1E6 # maximum possible data qubits (if number of data qubits is less than this, leave this as is)

qubitsorter = lambda q: (1 if q.qubit_type=='data' else MAXDATA) + q.qubit_id

def read(circ, with_wait=True):
	global MAXDATA
	def packable(qubits, used, total): # qubits = current qubits, used = current used qubits
		need = total[total.index(qubits[0]): total.index(qubits[-1])+1]
		return True not in [u in need for u in used] 

	def fillused(qubits, used, total):
		start = total.index(qubits[0])
		end = total.index(qubits[-1])
		return set(list(used) + total[start:end+1])

	qids = sorted([q.qubit_id for q in list(circ.qubit_gates_map.keys())])
	qubits = sorted(list(circ.qubit_gates_map.keys()), key=qubitsorter)

	ggrid = [circ.qubit_gates_map[q][:] for q in qubits]
	gnow = [ggrid[i].pop(0) if len(ggrid[i])>0 else [] for i in range(len(qids)) ]	

	while any(gnow):
		current = [None for i in qids]
		used = []
		for i,g in enumerate(gnow):
			if not g:  # no more gate on this qubit
				continue
			current_qubits = sorted(g.qubits[:], key=qubitsorter)
			if len(g.qubits) > gnow.count(g):	# not all ready
				continue
			if not packable(current_qubits, used, qubits): 
				continue
			indeces = [j for j in range(len(gnow)) if gnow[j] == g]
			current[i] = g	
			used = fillused(current_qubits, used, qubits)
			for j in indeces:
				gnow[j] = ggrid[j].pop(0) if ggrid[j] else None
		if not current:
			return
		yield current


def compute_locations(qubits, time):
	# qubits = [], time = int
	global qid_to_y, time_to_x, qsize, tsize
	qid_to_y = {}
	time_to_x = {}
	for i,q in enumerate(qubits):
		qid_to_y[q] = boxy*(2*i+1) + separationy*i
	for t in range(time+1):
		time_to_x[t] = boxx*(2*t+1) + separationx*(t+1)
	qsize = qubits
	tsize = time

def file2circuit(filename, with_latency=True):
	circ = circuit.Circuit()
	circ.load_from_file(filename, with_latency)
	return circ 

def draw_lines(qubits, total_time):
	startx = 0
	length = time_to_x[total_time]
	s = '\n'.join([add_line(startx, qid_to_y[q], length) for q in qubits])
	return s

def classical_line(startx, starty, length=10):
	style = styler(startx, starty-linewidth*2, length, linewidth)
	s = '<div class="qwire" %s>'%style + '</div>'
	style = styler(startx, starty-linewidth, length, linewidth*2)
	s += '<div class="qwire white_wire" %s>'%style + '</div>'
	style = styler(startx, starty+linewidth*2, length, linewidth)
	s += '<div class="qwire" %s>'%style + '</div>'
	return s

def draw_gate(gate, time, total_qubits, color='#FFF'):
	global linkmax

	qubits = [q for q in gate.qubits]
	x = time_to_x[time]
	y = qid_to_y[qubits[0]]
	boxsize = int(boxx*2)
	s = ''
	gatename = gate.gate_name.upper()

	if hasattr(gate, 'circuit_list'):
		link = '%d.html'%linkmax
		linkmax += 1
		links[gate] = link
	else:
		link = ''			

	if len(qubits)==1:
		if 'M' == gatename or 'METER' in gatename or 'MEASURE' in gatename:
			s += classical_line(x, y, time_to_x[tsize]-x)		
		if 'CORRECT' in gatename: 
			gatename = 'C'	
		if hasattr(gate, 'is_error') and gate.is_error:
			color = '#F99'
		s += gate_drawer.draw_one_q_gate(gate, x, y, boxsize, boxsize, color, link)
	else:
		ys = [qid_to_y[q] for q in qubits]
		unused_ys = [qid_to_y[q] for q in total_qubits 
			if (q not in qubits and qid_to_y[q] > min(ys) and qid_to_y[q] < max(ys))]
		s += gate_drawer.draw_multi_q_gate(gate, x, ys, unused_ys, boxsize, boxsize, color, link)
	return s

def draw_qubits(qubits):
	qubits = sorted(qubits)
	s = '\n'.join([add_label('%s %d'%(q.qubit_type[0], q.qubit_id), 0, qid_to_y[q], boxx*2, boxy*1) for q in qubits])
	return s

def circ_str(circ, qubits, draw_also=[]):
	s = ''
	circ_array = list(read(circ))
	compute_locations(qubits, len(circ_array))
	s += draw_lines(qubits, len(circ_array))
	t = 0
	for gates in circ_array:
		for gate in gates:
			if gate:
				s+= '%s \n'%draw_gate(gate, t, qubits)
		s += '<br>'
		#t -= 1
		t += 1
	return s

def draw_circuit(circ, offset=0, draw_also=[]):
	qubits = sorted([q for q in list(circ.qubit_gates_map.keys()) if not q.level], key=qubitsorter)
	circstr = circ_str(circ, qubits, draw_also)
	qubitstr = draw_qubits(qubits)
	x = max(time_to_x.values())
	y = max(qid_to_y.values())	
	qubit_area_width = 60
	style = styler(5, 45, qubit_area_width, y+boxy*10)
	s = '<div class="qubitarea" %s>'%style
	s += qubitstr
	s += '</div>'
	style = 'style="position: absolute; left: %dpx; width:%dpx; height:%dpx;"'%(
			qubit_area_width, x+boxx*len(qubits), y+boxy*10)
	s += '<div class="circuitarea" %s>'%style
	s += circstr
	s += '</div>'
	return s

def draw_circuits(circuit_list, draw_also=[]):
	global time_to_x
	s = ''
	for circ in circuit_list:
		# XXX edit here for use actual time
		# if 'time' in draw_also: blah
		offset = max(time_to_x.values())+20 if time_to_x else 0 
		s += draw_circuit(circ, offset, draw_also=[])
	return s

def draw_circuits_recursively(circs, parent_gate_name='main', link='main.html', draw_also=[]):
	circ_html = draw_circuits(circs)
	yield circ_html, parent_gate_name, link
	for circ in circs:
		for g in circ.gates:
			if g.qubits and hasattr(g, 'circuit_list'):
				if not g.circuit_list:
					link = links[g]
					yield 'empty', 'this is empty', link
					continue
				if len(g.circuit_list) > 1:
					print(' * New! * many circuits in a list in ', g.gate_name)
					print(' * Check if offset feature is working well *')
				for child in draw_circuits_recursively(g.circuit_list, g.gate_name, links[g]):
					yield child
	
def draw_from_file(filename=None):
	if not filename.strip():
		filename = 'steaneT.txt'
	circ = file2circuit(filename)
	return draw_from_circuit(circ)

