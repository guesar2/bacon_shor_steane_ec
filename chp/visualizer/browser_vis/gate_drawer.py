from .htmler import *
import os
linewidth = 1

path = os.path.dirname(__file__)
path = '../'

logger = 'log.txt'

gate_images = {'M': path+'images/meter.png'}

def remove_line(startx, starty, length=10):
	return add_line(startx, starty, length, linewidth, 'white_wire')

def draw_a_box(x, y, sizex, sizey, label='', cls='', extra=''):
	style = styler(x-sizex//2, y-sizey//2, sizex, sizey, extra)
	return '<div class="gate %s" %s><p class="gate">%s</p></div>'%(cls, style, label)

def write_annotation(s, x, y, sizex, sizey, vertical=False, link=''):
	if vertical:
		cls = 'verticallabel'
		style = styler(x-sizex//2, y, sizex, '')	
	else:
		cls = 'annotation'
		style = styler(x-sizex//2, y, sizex, sizey)
	st = '<a href="%s"> %s </a>'%(link, s) if link else ' %s '%s
	return 	'<p class="%s" %s> %s </p>'%(cls, style, st)

def draw_target(x, y, size):
	size = size*3//4
	s = add_circle(x-size//2+1, y-size//2, size, 'target')
	s += add_line(x+1, y-size//2, linewidth, size+1, 'vertical_wire')
	return s

def draw_control(x, y, size, color='#000', extra=''):
	s = add_circle(x-size//2, y-size//2, size, 'control')
	return s

def draw_container_gate(name, x, ys, nys, xsize, ysize, link=''):
	xsize = xsize*3//4
	ysize = ysize*3//4	
	name = name.upper()
	ysi = min(ys)
	ysf = max(ys)
	dys = ysf - ysi
	s = add_line(x+1, ysi, linewidth, dys)
	s+= draw_a_box(x, ysi+dys//2, xsize, dys+ysize, '', 'container')#name[0])
	s+= write_annotation(name, x, min(ys), xsize, ysize, True, link)
	radius = 2
	for y in ys: # used qubits y values
		s += add_circle(x-xsize//2-radius, y-radius, radius*2, 'container_connection')
		s += add_circle(x+xsize//2-radius+2, y-radius, radius*2, 'container_connection')
	for y in nys: # unused qubits y values
		s += add_line(x-xsize//2, y, xsize, 1, 'dashed_wire')
	return s

def draw_multi_q_gate(gate, x, ys, unused_ys, xsize, ysize, color='#000', link=''):
	name = gate.gate_name.upper()
	if 'C' == name:
		name = 'CX'
	if 'CX' in name and name.strip('X') == 'C':
		name = 'CX'
	if 'CZ' in name and name.strip('Z') == 'C':
		name = 'CZ'
	if not (name == 'CX' or name == 'CZ'): # not CX nor CZ -> container gate
		return draw_container_gate(name, x, ys, unused_ys, xsize, ysize, link)
	xsize = xsize*3//4
	ysize = ysize*3//4	
	name = name.upper()
	s = add_line(x+1, min(ys), linewidth, max(ys)-min(ys), 'vertical_wire')
	
	yi = min(ys)-ysize*0.65
	yf = max(ys)+ysize*0.48
	s += get_time_strings(gate, x, yi, yf, xsize, ysize)
	
	for i, y in enumerate(ys):
		if i == 0:
			s += draw_control(x, ys[i], 8, color)
		elif name[1] == 'X':
			s+= draw_target(x, ys[i], ysize)
		else:
			s+= draw_a_box(x, ys[i], xsize, ysize, name[1])
	return s

def get_time_strings(gate, x, yi, yf, xsize, ysize):
	s = ''
	if hasattr(gate, 'start_time') and gate.start_time != None:		
		s += write_annotation(str(gate.start_time), x, yi, xsize, ysize)
		s += write_annotation(str(gate.end_time), x, yf, xsize, ysize)
	return s

def draw_one_q_gate(gate, x, y, xsize, ysize, color='#FFF', link=''):
	color = "background-color: %s;"%color

	gatename = gate.gate_name.upper()
	s = ''
	yi = y-ysize*0.65
	yf = y+ysize*0.48
	s += get_time_strings(gate, x, yi, yf, xsize*3//4, ysize*3//4)	

	if 'M' == gatename or 'METER' in gatename or 'MEASURE' in gatename:
		addr = gate_images['M']
		s += add_pict(addr, x-xsize//2, y-ysize//2, xsize, ysize)
		return s

	ysize = ysize*3//4

	xsize = xsize*3//4

	if 'PREPARE' in gatename:
		label = 'prep'
		nowirefrom = x-xsize*3//2
		s += remove_line(nowirefrom, y, x + xsize - nowirefrom)
		s += draw_a_box(x, y-ysize*0.3, xsize, ysize//2, label, 'prepare')
		s += add_pict(path+'images/state_%s.png'%gatename[7:], x-xsize, y+ysize*0.1, xsize*2, ysize//2)
		return s
	if 'WAIT' in gatename:
		label = 'I' # gate.name[0]
		s += draw_a_box(x, y, xsize, ysize, label, 'light')
		s += write_annotation(name[4:], x, y, xsize, ysize)
		return s

	return s+draw_a_box(x, y, xsize, ysize, gatename[0].upper(), '', color)
