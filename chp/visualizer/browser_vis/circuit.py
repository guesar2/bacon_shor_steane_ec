#!/usr/bin/python

# circuit.py
# This is for tempolary use and only used when chp files are read by local.py
# This file will be removed once load_from_file() is integrated into local or circuit_drawer.
# -Yu

class Qubit(object):
	def __init__(self, qid):
		self.qubit_id = qid

class Gate(object):
	def __init__(self, name, qubits):
		self.gate_name = name
		self.qubits = qubits
	def __repr__(self):
		return 'GATE: %s, %r'%(self.gate_name, self.qubits)

class Circuit(object):
	def __init__(self, gates=[], qids=[]):
		self.gates = gates
		self.qubit_gates_map = None
		self.update_gate_map()

	def update_gate_map(self):
		m = {}
		for g in self.gates:
			for q in g.qubits:
				m[q] = m.get(q, []) + [g]
		self.qubit_gates_map = m

	def load_from_file(self, filename, with_latency=True):
		translater = {'C': 'CX', 'CNOT': 'CX'}
		f = open(filename)
		created_qubits = {}
		for line in f:
			line = line.replace(',',' ').split()
			if not line:
				continue
			if '#' in line[0]:
				continue
			name = line[0].upper()
			name = translater.get(name, name)
			if '.' in line[-1]:
				latency = line.pop()
				if with_latency:
					name += ('('+latency+')')
			qids = [int(q) for q in line[1:]]
			qubits = []
			for q in qids:
				try:
					qubits.append(created_qubits[q])
				except KeyError:
					new_qubit = Qubit(q)
					qubits.append(new_qubit)
					created_qubits[q] = new_qubit
			self.gates.append(Gate(name, qubits))
		self.update_gate_map()
