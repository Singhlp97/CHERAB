loading=['|', '/', '-', '\\']

def progress_bar_Ricky(i,num_tot):
	pb_len = 50
	step=int(round((i+1)/num_tot*pb_len))
	if (i+1)<num_tot:
		print('\r     Progress: [' + '#'*step + ' '*(pb_len-step) + '] ' + '{} %'.format(round((i+1)/num_tot*100,1)) + ' {}'.format(loading[i%4]), end='\r')
	else:
		print('\r     Progress: [' + '#'*step + ' '*(pb_len-step) + '] ' + '{} %'.format(round((i+1)/num_tot*100,1)) + ' Done!', end='\n\n')

