from main import *


# In[5]:


identification='slopes'
# identification='height'
if identification=='slopes':
    power=2
elif identification=='height':
    power=0


# In[7]:



if __name__ == "__main__":
    KT=array([0.02,2000])
    U10=5
      
    x=linspace(0,400,200)
    y=linspace(0,400,200)
    k=logspace(log10(KT[0]),log10(KT[-1]),128)
    phi=linspace(-pi,pi,50)
    t=linspace(0,120,1)

    m=model(k,phi,t,U10,KT)
    model_t = m.main(k,phi,t) 

    x, y = np.meshgrid(x, y)
    
    
    for i in range(len(t)):
        figure()
        z=model_t[i]([x,y],phi,t=[t[i]])
        contourf(z,100,cmap=cm.winter)
        colorbar()
        ylabel(r'Y, \text{м}',fontsize=16)
        xlabel(r'X, \text{м}',fontsize=16)
    #   savefig(path.abspath('..'+'\\water\\anim\\'+'water'+str(i)+'.png'),
    #             pdi=10**6,bbox_inches='tight')
        show()

