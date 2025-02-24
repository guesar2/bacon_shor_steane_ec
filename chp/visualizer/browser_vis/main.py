#!/usr/bin/python

import cgi
import cgitb
from . import circuit_drawer
import os
import json
cgitb.enable()

def typehtml():
	print("Content-type: text/html\n")

def typejson():
	print('Content-type: application/json\n')	

def params():
	files = [f for f in os.listdir('inputs') if f[0] != '.']
	return json.dumps({'filenames':files})

def main():
	form = cgi.FieldStorage()	

	if 'param_request' in form:
		typejson()
		print(params())
		return

	typehtml()

	if 'filename' in form:
		filename = form.getvalue('filename')
		if filename == 'undefined':
			filename = None
		print(circuit_drawer.draw(filename))

main()
