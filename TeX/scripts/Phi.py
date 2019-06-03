from main import *
rcParams['figure.figsize'] = [4, 4]
rcParams['axes.labelsize'] = 20
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


k_m=m.k_m
x=linspace(-pi,pi,10000)
temp=[k_m/2,k_m,2*k_m]
for t in temp:
    y=m.Phi(t,x)
    polar(x,y/k_m,label=r'$k/k_m=$' + str(round(t/k_m,3)))
legend()
savefig(path.abspath('..'+'\\fig\\full_angles1.pdf'),
       bbox_inches='tight')
show()
temp=[10*k_m,50*k_m,100*k_m]
for t in temp:
    y=m.Phi(t,x)
    polar(x,y/k_m,label=r'$k/k_m=$' + str(round(t/k_m,3)))
legend()
savefig(path.abspath('..'+'\\fig\\full_angles2.pdf'),
       bbox_inches='tight')
show()