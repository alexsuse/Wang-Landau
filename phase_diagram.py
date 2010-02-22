import os
import multiprocessing
from multiprocessing import Process, Queue
from Queue import Empty

def runWS(q):
	while True:
		try:
			x = q.get(block=False)
			
			[RP,D,L,flatness,fmin,P,delta,thread,irun] = x
			nome='_d'+str(int(delta*10))+'_P'+str(P)+'_L'+str(L)+'_D'+str(D)+'_run'+str(irun)+'.dat '
	  
		      	os.system('./moralWLwindows '+str(delta)+' '+str(P)+' '+str(RP)+' '
			+str(D)+' '+str(L)+' '+str(flatness)+' '+str(fmin)
			+' '+'wlsrange'+str(thread)+nome+'dos'+str(thread)+nome+'thermo'
			+str(thread)+nome+'canonica'+str(thread)+nome+'>log'+str(thread)+nome)
			
			print 'thread', thread,' run', irun, ' P=', P, ' ', ' delta=', delta				
			
		except Empty:
			break	


if __name__ == '__main__':

	#cria nthreads (=numero de cores)
	nthreads=multiprocessing.cpu_count()
	
	#parametros fixos
	D=5  		 # dimensao da matriz moral
	RP=0		 # rewiring probability
	L=10		 # tamanho da rede quadrada LxL
	flatness=0.85    # flat histogram flatness
	fmin=0.00001	 # minimum f 
	par_fixos=[RP,D,L,flatness,fmin]

	#prepara fila de tarefas
	work_queue = Queue()
	for copia in range(0,1,1):
		for irun in range(0,5,1):
			for P in range(1,5,1): 
				delta=1.0
				while(delta>0.05):			
					work_queue.put(par_fixos+[P,delta,copia,irun])
					delta=delta-0.1

	#multiprocessamento
	processes = [Process(target=runWS, args=(work_queue,)) for i in range(nthreads)]		
	
	for p in processes:
		p.start()
	for p in processes:
		p.join()

	#cria diretorio para dados	
	lista= os.listdir(os.getcwd())
	if not ('log' in lista):os.system('mkdir logs')
	if not ('WLSrange' in lista):os.system('mkdir WLSrange')
	if not ('thermo' in lista):os.system('mkdir thermo')
	if not ('canonical' in lista):os.system('mkdir canonical')
	if not ('dos' in lista):os.system('mkdir dos')

	#transfere resultados
	os.system('mv log*.dat log')
	os.system('mv wlsrange*.dat log')
	os.system('mv thermo*.dat log')
	os.system('mv canonica*.dat log')
	os.system('mv dos*.dat log')


