from main import *
rcParams['figure.figsize'] = [6, 6]
rcParams['axes.labelsize'] = 18
rc('text', usetex=True)
rc('text.latex', preamble=[r'\usepackage[russian]{babel}',
                         r'\usepackage{amsmath}',
                        r'\usepackage{amssymb}',
                        r'\usepackage{mathrsfs}'])
rc('font', family='serif')
if __name__ == "__main__":
    KT=array([0.02,2000])
    U10=5
      
    x=linspace(0,400,200)
    y=linspace(0,400,200)
    k=logspace(log10(KT[0]),log10(KT[-1]),128)
    phi=linspace(-pi,pi,50)
    t=linspace(0,120,1)

    m=model(k,phi,t,U10,KT)

x=20170

U10=5
k=logspace(log10(KT[0]),log10(KT[-1]),1000)
C=k**2*m.spec.full_spectrum(k,x)
loglog(k,C,label=r'$U_{10}=5 \text{ м/с}$')

U10=10
m=model(k,phi,t,U10,KT)
k=logspace(log10(KT[0]/10),log10(KT[-1]),1000)
C=k**2*m.spec.full_spectrum(k,x)
loglog(k,C,label=r'$U_{10}=10 \text{ м/с}$')

U10=15
m=model(k,phi,t,U10,KT)
k=logspace(log10(KT[0]/10),log10(KT[-1]),1000)
C=k**2*m.spec.full_spectrum(k,x)
loglog(k,C,label=r'$U_{10}=15 \text{ м/с}$')

U10=20
m=model(k,phi,t,U10,KT)
k=logspace(log10(KT[0]/10),log10(KT[-1]),1000)
C=k**2*m.spec.full_spectrum(k,x)
loglog(k,C,label=r'$U_{10}=20 \text{ м/с}$')
legend(fontsize=16)

xlabel(r'$k, \text{рад}\cdot\text{м}^{-1}$')
ylabel(r'$S_{\theta}$')
ylim((10**(-15),1000))
savefig(path.abspath('..'+'\\fig\\full_spectrum3.pdf'),
       bbox_inches='tight')
show()


xx=arange(1430,20000,4000)
xx=hstack((xx,[20170]))
U10=10
for n in xx:
    m=model(k,phi,t,U10,KT,n)
    k=logspace(log10(k[0]),log10(k[-1]),1000)
    C=k**2*m.spec.full_spectrum(k,n)
    loglog(k,C,label=r'$\widetilde{x}=$'+'{0}'.format(n))

# ylim([10**-30,10])
xlabel(r'$k, \text{рад}\cdot\text{м}^{-1}$')
ylabel(r'$S_{\theta}$')
legend(fontsize=16)
ylim((10**(-15),1000))
savefig(path.abspath('..'+'\\fig\\full_spectrum4.pdf'),
       bbox_inches='tight')

show()