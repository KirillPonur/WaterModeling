from main import *
rcParams['figure.figsize'] = [6, 4]
rcParams['axes.labelsize'] = 14
rc('text', usetex=True)
rc('text.latex', preamble=[r'\usepackage[russian]{babel}',
                         r'\usepackage{amsmath}',
                        r'\usepackage{amssymb}',
                        r'\usepackage{mathrsfs}'])
rc('font', family='serif')


if __name__ == "__main__":

    identification='slopes'
    # identification='height'
    if identification=='slopes':
        power=2
    elif identification=='height':
        power=0


    KT=array([0.02,2000])
    
    x=linspace(0,400,200)
    y=linspace(0,400,200)
    k=logspace(log10(KT[0]),log10(KT[-1]),128)
    phi=linspace(-pi,pi,50)
    t=linspace(0,120,1)
    rho=linspace(0,100,100)
    N=25
    M=50
    U10=8

    m=main(N,M,t, U10, KT)
    model_t = m.main(m.k,m.phi,t)

      






ki,b0,err=m.interspace(m.k0,int(N),cache=False)
ki=m.nodes(ki,b0,cache=False)


# # In[38]:

y=m.full_spectrum(m.k0)
for i in range(len(ki)):
    loglog([ki[i],ki[i]],[KT[0]**2*m.full_spectrum(KT[0]),ki[i]**2*m.full_spectrum(ki[i])],color='darkblue')    
loglog(m.k0,m.k0**2*y,color='crimson')
# legend(fontsize=16)
xlabel(r'$k, \text{рад} \cdot\text{м}^{-1}$',fontsize=16)
ylabel(r'$S_{\theta},\text{a.u.}$',fontsize=16)
grid(which='major', linestyle='-')
grid(which='minor', linestyle=':')
minorticks_on()

savefig(path.abspath('..'+'\\fig\\split_angles'+'.pdf'),
       bbox_inches='tight')
show()


