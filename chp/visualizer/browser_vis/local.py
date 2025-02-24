from . import circuit_drawer
import subprocess
import os
import sys
import platform

path = os.path.abspath(__file__)
path = path.strip('local.py')
directory = os.path.dirname(__file__)
tempdir = directory + '/temp/'

def header(title=''):
	s = '''<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">

<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
	<title>%s</title>
	<script type="text/javascript" src="jquery-1.7.1.min.js"></script>	
	<link rel="stylesheet" type="text/css" href="style.css" media="all" />	
  <body>
	<div id="main">
'''%(title)
	return s

def footer(comment):
	s = '''	</div>
	<div id="header">
	%s
	</div>			
	<div id="footer">
		<form method="post">
			<div id='form_elements'>
			</div>
		</form>
	</div>	
  </body>
  </html>'''%comment
	return s


def make_page(content, location='tempcircuit.html', title="", comment=""):
	filename = tempdir+location
	f = open(filename, 'w')
	f.write(header(title))
	f.write(content)
	f.write(footer(comment))
	f.close()
	return filename

def show(filename):
	mac = ('Darwin' == platform.system())
	if not mac:
		try:
			subprocess.Popen(['google-chrome', '%s'%(filename)])
		except OSError:
			print('* google-chrome not found. Trying with firefox...')
			try:
				subprocess.Popen(['firefox', '%s'%(filename)])
			except OSError:
				print('* firefox not found. Please try again after installing firefox or google-chrome. Aborting.')
				sys.exit(1)
	else:
		subprocess.Popen(['open', '%s'%filename])

def from_file(filename):
	'''One of main public functions. Use this to visualize a circuit from file.
	Use from_circuit(circ) instead to visualizer from a python circuit object.'''
	f = open(tempfile, 'w')
	f.write(header())
	f.write(circuit_drawer.draw_from_file(filename))
	f.write(footer())
	f.close()
	os.system('google-chrome %s'%tempfile)

def from_circuit(circ, recursive=False, show_error=False, show_time=False):
	'''One of main public functions. Use this to visualize from a python circuit object.
	Use from_file(filename) instead to visualize from a text file with list of gates'''
	draw_also = []
	if show_error:
		draw_also += ['error']
	if show_time:
		draw_also += ['time']
	if recursive:
		for s, title, link in circuit_drawer.draw_circuits_recursively([circ], draw_also=draw_also):
			make_page(s, link, title, title)	
			if title == 'main':
				show(tempdir+link)
	else:
		s = circuit_drawer.draw_circuit(circ, draw_also=draw_also)
		filename = make_page(s)
		show(filename)

def debug(circ):
	qubitsorter = lambda q: (1 if q.qubit_type=='data' else 1E6) + q.qubit_id
	qubits = sorted(list(circ.qubit_gates_map.keys()), key=qubitsorter)
	print(qubits)
