from Tkinter import *
import os

root = Tk()               # root (main) window
top = Frame(root)         # create frame
top.pack(side='top')      # pack frame in main window

hwtext = Label(top, text='Wang Landau    ')
hwtext.pack(side='left')

hwtext = Label(top, text='delta')
hwtext.pack(side='left')
delta = StringVar() # variable to be attached to delta_entry
delta.set('1.0')    # default value
delta_entry = Entry(top, width=4, textvariable=delta)
delta_entry.pack(side='left')

hwtext = Label(top, text='P')
hwtext.pack(side='left')
P = StringVar() # variable to be attached to P_entry
P.set('1')    # default value
P_entry = Entry(top, width=4, textvariable=P)
P_entry.pack(side='left')

hwtext = Label(top, text='L')
hwtext.pack(side='left')
L = StringVar() # variable to be attached to L_entry
L.set('10')    # default value
L_entry = Entry(top, width=4, textvariable=L)
L_entry.pack(side='left')

hwtext = Label(top, text='D')
hwtext.pack(side='left')
D = StringVar() # variable to be attached to D_entry
D.set('5')    # default value
D_entry = Entry(top, width=4, textvariable=D)
D_entry.pack(side='left')

hwtext = Label(top, text='RP')
hwtext.pack(side='left')
RP = StringVar()
RP.set('0.0')    
RP_entry = Entry(top, width=4, textvariable=RP)
RP_entry.pack(side='left')

hwtext = Label(top, text='flatness')
hwtext.pack(side='left')
flatness = StringVar() 
flatness.set('0.85')  
flatness_entry = Entry(top, width=4, textvariable=flatness)
flatness_entry.pack(side='left')

hwtext = Label(top, text='minlogf')
hwtext.pack(side='left')
minlogf = StringVar() 
minlogf.set('0.0001')  
minlogf_entry = Entry(top, width=10, textvariable=minlogf)
minlogf_entry.pack(side='left')

hwtext = Label(top, text='irun')
hwtext.pack(side='left')
irun = StringVar() 
irun.set('0')  
irun_entry = Entry(top, width=4, textvariable=irun)
irun_entry.pack(side='left')

hwtext = Label(top, text='thread')
hwtext.pack(side='left')
t = StringVar() 
t.set('0')  
t_entry = Entry(top, width=4, textvariable=t)
t_entry.pack(side='left')

def runWL():
    deltax=float(delta.get()) 	
    Px=int(P.get())
    Lx=int(L.get())
    Dx=int(D.get())
    irunx=int(irun.get())
    RPx=float(RP.get())
    flatnessx=float(flatness.get())
    fminx=float(minlogf.get())
    threadx=int(t.get())  	
	
    nome='_d'+str(int(deltax*10))+'_P'+str(Px)+'_L'+str(Lx)+'_D'+str(Dx)+'_run'+str(irunx)+'.dat '
  
    os.system('./moralWLwindows '+str(deltax)+' '+str(Px)+' '+str(RPx)+' '
    +str(Dx)+' '+str(Lx)+' '+str(flatnessx)+' '+str(fminx)
    +' '+'wlsrange'+str(threadx)+nome+'dos'+str(threadx)+nome+'thermo'
    +str(threadx)+nome+'canonica'+str(threadx)+nome)

    print 'thread', t,' run', irun, ' P=', P, ' ', ' delta=', delta	

compute = Button(top, text=' run  ', command=runWL)
compute.pack(side='left')
root.mainloop()

