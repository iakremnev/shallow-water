import math

x_cells = 100
y_cells = 100

prefix = 'input/'

def gen_height():

	f = open(prefix + 'height.txt', 'w')	
	f.write(str(x_cells) + ' ' + str(y_cells) + '\n')
	for i in range(x_cells):
		for j in range(y_cells):
			#f.write('1.5 ')
			f.write('3.0 ' if (i < 20 and j < 20) else '1.0 ')
			#f.write('2.0 ' if (i*i*0.25 < j) else '1.0 ')
			#f.write(str(math.exp(-((i-x_cells/2)**2 + (j - y_cells/2)**2)/100) + 1.0) + ' ')
			#f.write('3.0 ' if ((i - j)**2 <= 16) else (str(abs(math.sqrt(i*i + 100)/100 - math.sqrt(j))) + ' '));
		f.write('\n')
		
			
def gen_velocity_x():

	f = open(prefix + 'velocity_x.txt', 'w')
	f.write(str(x_cells + 1) + ' ' + str(y_cells) + '\n')
	for i in range(x_cells + 1):
		for j in range(y_cells):
			f.write('0.0 ')
		f.write('\n')

		
def gen_velocity_y():

	f = open(prefix + 'velocity_y.txt', 'w')
	f.write(str(x_cells) + ' ' + str(y_cells + 1) + '\n')
	for i in range(x_cells):
		for j in range(y_cells + 1):		
			f.write('0.0 ')
		f.write('\n')
		
		
def gen_surface():

	f = open(prefix + 'surface.txt', 'w')
	f.write(str(x_cells) + ' ' + str(y_cells) + '\n')
	for i in range(x_cells):
		for j in range(y_cells):
			#f.write('0.0 ')
			f.write(str(math.exp(-((i-x_cells/2)**2 + (j-y_cells/2)**2)/100)) + ' ')
		f.write('\n')

gen_height()
gen_velocity_x()
gen_velocity_y()
gen_surface()
