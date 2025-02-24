def styler(x, y, width='', height='', extra=''):
	if not width:
		width = 'auto'
	else:
		width = str(width) + 'px'
	if not height:
		height = 'auto'
	else:
		height = str(height) + 'px'
	style = 'style="width:%s; height:%s;'%(width, height)
	style+= 'position:absolute; top:%dpx; left:%dpx;'%(y, x)
	style+= '%s"'%extra
	return style

def add_pict(url, x, y, width, height, extra=''):
	style = styler(x, y, width, height, extra)
	return '<img src="%s" %s />'%(url, style)

def add_circle(x, y, size, cls=''):
	extra = '-moz-border-radius: %dpx; -webkit-border-radius: %dpx;'%(size, size)
	style = styler(x, y, size, size, extra)
	return '<div class="circle %s" %s > </div>'%(cls,style)

def add_line(x, y, length=100, thick=1, cls=''):#color='#000', decor='solid'):
	style = styler(x, y, length, thick)
	return '<div class="qwire %s" %s>'%(cls, style) + '</div>'

def rotation_style():
	style = '-moz-transform: rotate(90deg);' # firefox
	style +='-webkit-transform: rotate(90deg);' 
	style += 'filter: progid:DXImageTransform.Microsoft.BasicImage(rotation=1);' # IE
	return style

def add_label(s, x, y, sizex, sizey):
	style = styler(x, y, sizex, sizey)
	return '<p %s>%s</p>'%(style, s)
