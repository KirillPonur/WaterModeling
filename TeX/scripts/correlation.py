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
    N=256
    M=50
    U10=8

    m=main(N,M,t, U10, KT)
    model_t = m.main(m.k,m.phi,t)

      

ki,b0,err=m.interspace(m.k0,int(2*N),cache=False)
ki=m.nodes(ki,b0,cache=False)

# In[29]:
rho=linspace(0,400,400)
y=zeros(len(rho))
k=logspace(log10(m.KT[0]),log10(m.KT[-1]),10**5)
y1=m.height(k,rho)
k=logspace(log10(m.KT[0]),log10(m.KT[-1]),256)
# k=linspace(m.KT[0],m.KT[-1],256)
y=m.height_sum(k,rho)
y2=m.height_sum(ki,rho)
y2=y2/max(y2)
y1=y1(rho)
real=y1/max(y1)
modeling=y/max(y)
plot(rho,modeling,label='Логарифмическое поле',color='darkgreen')
plot(rho,y2, color='blue',label='<<Отбеленное>> поле')
plot(rho,real,label='Реальное поле',color='crimson')
legend(fontsize=16)
xlabel(r'${\rho},\text{м}$',fontsize=16)
ylabel(r'$K, a.u.$',fontsize=16)
grid(which='major', linestyle='-')
grid(which='minor', linestyle=':')
minorticks_on()

savefig(path.abspath('..'+'\\fig\\correlation_height_wa'+'.pdf'),
       bbox_inches='tight')
show()






# In[38]:


# rho=linspace(0,400,400)
# y=zeros(len(rho))
# k=logspace(log10(m.KT[0]),log10(m.KT[-1]),10**5)
# y1=m.height(k,rho)
# y=m.height_sum(ki,rho)
# y1=y1(rho)
# real=y1/max(y1)
# modeling=y/max(y)
# plot(rho,modeling,label='Модельное поле',color='blue')
# plot(rho,real,label='Реальное поле',color='crimson')
# legend(fontsize=16)
# xlabel(r'${\rho},\text{м}$',fontsize=16)
# ylabel(r'$K_{\theta}, a.u.$',fontsize=16)
# grid(which='major', linestyle='-')
# grid(which='minor', linestyle=':')
# minorticks_on()

# # savefig(path.abspath('..'+'\\fig\\correlation_height_wa'+'.pdf'),
# #        bbox_inches='tight')
# show()


rho=linspace(0,400,400)
y=zeros(len(rho))
k=logspace(log10(m.KT[0]),log10(m.KT[-1]),10**5)
y1=m.angles(k,rho)
k=logspace(log10(m.KT[0]),log10(m.KT[-1]),256)
# k=linspace(m.KT[0],m.KT[-1],256)
y=m.angles_sum(k,rho)
y2=m.angles_sum(ki,rho)
y2=y2/max(y2)
y1=y1(rho)
real=y1/max(y1)
modeling=y/max(y)
plot(rho,modeling,label='Логарифмическое поле',color='darkgreen')
plot(rho,y2,color='blue',label='<<Отбеленное>> поле')
plot(rho,real,label='Реальное поле',color='crimson')
legend(fontsize=16)
xlabel(r'${\rho},\text{м}$',fontsize=16)
ylabel(r'$K, a.u.$',fontsize=16)
grid(which='major', linestyle='-')
grid(which='minor', linestyle=':')
minorticks_on()

savefig(path.abspath('..'+'\\fig\\correlation_angles_wa'+'.pdf'),
       bbox_inches='tight')
show()


