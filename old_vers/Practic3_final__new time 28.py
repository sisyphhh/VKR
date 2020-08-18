#Allel`s ficsation in 9 sysytem: 0 or 10 "mutation" fixsation
import numpy as np
from matplotlib import pylab as plt
from matplotlib import rcParams
get_ipython().run_line_magic('matplotlib', 'inline')
import random
import os
import pandas as pd
from scipy.stats import gamma
from scipy.stats import erlang
# ### Условия на лямбду(на матрицу переходов)
#     λ(0,1)=2N*μ(2,1) λ(0,0)=-2N*μ(2,1) 
#     λ(N,N-1)=2N*μ(1,2) λ(N,N)=-2N*μ(1,2) 
#     λ(i,i-1)=[i(2N-i)/(2N)] + μ(2,1)*(2N-i)     
#     λ(i,i+1)=[i(2N-i)/(2N)] + μ(1,2)*(i)        
#     λ(i,i)=-λ(i,i-1)-λ(i,i+1)

#необходимые переменные; общий случай
eps=0.00001
NNN=5
liamda_i_massive=[]
matrix_liambda=[]
def catalogmaker(name='Graphics'):
    if not os.path.exists(name): 
        os.mkdir(name)# создание каталога 
        print('Каталог успешно создан', name)
    main_target=os.getcwd() + os.sep + name + os.sep
    return main_target

def matr_printer(matrix_liambda,mu_12,mu_21=-1,NNN=5):#просто печать/вывод матрицы(списка значений)
    if mu_21==-1:
        print('\nLiambda matrix A_μ(1,2)==μ(2,1)=="%s":'%mu_12) #μ(1,2)≠μ(2,1)
    else:
        print('\nLiambda matrix A_μ(1,2)="%s" ;μ(2,1)="%s":'%(mu_12,mu_21))
    for i in range(2*NNN +1):
        print()
        for j in range(2*NNN +1):
            print(matrix_liambda[i][j],' ',end='')
    print('\n\n')

#создание матрицы интенсивностей (нужно добавить доп параметры на mu_12 mu_21когда не раваны)
def creator_mart_liambd(mu_12,mu_21=-1,NNN=5):   
    if mu_21==-1:#на случай дефолта
        mu_21=mu_12
    
    matr_liam=[] 
    for i in range(2*NNN +1):
        matr_liam.append([])
        #print(matrix_liambda)
        for j in range(2*NNN +1):
            if (i == 0):
                if j==i:
                    matr_liam[i].append(-2*NNN*mu_21)#μ(2,1)
                elif j==i+1:
                    matr_liam[i].append(2*NNN*mu_21)
                else:
                    matr_liam[i].append(0.0)
            elif (i==10):
                if j==i:
                    matr_liam[i].append(-2*NNN*mu_12)#μ(1,2)
                elif j==i-1:
                    matr_liam[i].append(2*NNN*mu_12)
                else:
                    matr_liam[i].append(0.0)
            else:
                if j==i:
                    matr_liam[i].append(round(-(i*(2*NNN-i)/(2*NNN) + mu_21*(2*NNN-i))-(i*(2*NNN-i)/(2*NNN) + mu_12*(i)),4))
                elif (j==i-1):    
                    matr_liam[i].append(round(i*(2*NNN-i)/(2*NNN) + mu_21*(2*NNN-i),4))#
                elif (j==i+1):  
                    matr_liam[i].append(round(i*(2*NNN-i)/(2*NNN) + mu_12*(i),4))#
                else:
                    matr_liam[i].append(0.0)
    return matr_liam


#3
def Graph_Ri_i(r_arr, Mu,Mu2=0, N=11, label_name='Dependence`s Graph of Ri from i=0..10',add=''):
    rcParams['figure.dpi']=150
    fig, ax = plt.subplots(figsize=(10,6))
    main_target=catalogmaker(name='Graphics')
    x_i=[x for x in range(0,11)]
    print(x_i)
    if Mu2==0:#те для симметричного варианта - просто такой дефолт 
        for i in range(len(r_arr)):
            plt.plot(x_i,r_arr[i],'--',label=Mu[i])
    else:
        for i in range(len(r_arr)):
            plt.plot(x_i,r_arr2[i],'--',label='mu_12="%s"---mu_21="%s"'%(Mu[i],Mu[2]))

    #plt.plot(x_i,y_i,'--',color='black',label='h1=0,1')
    plt.title('%s'%(label_name))
    plt.ylabel('Ri');plt.xlabel('Xi')
    plt.legend()
    if Mu2==0:#те для симметричного варианта - просто такой дефолт 
        plt.yticks(np.arange(0, np.max(r_arr)+0.05, 0.02))
    else:
        plt.yticks(np.arange(0, np.max(r_arr2)+0.05, 0.05))  
    plt.xticks(np.arange(0, 10+1, 1.0))
    plt.grid(True)
    #plt.show()
    fig.savefig(main_target+'stationar_'+add)

#r_arr=[[],[],...] 11 шт по 11 эл-тов; Graph_Ri_i(r_arr,Mu,add='simmetr')

# 4-5 
#функция считающая мат_ожидание особей по частотам, получ. из стационарных ур.
def average_funck(r_ar,x_i,mu,par='A1',NNN=5,len_ar=2):
    sum_C_R=[]
    if len_ar==1:
        box=0
        for i in range(len(r_ar)):
            if par=='A1':
                box+=x_i[i]*r_ar[i]
                #print(' ',x_i[j],'*',r_arr2[i][j],end='')
            else: #par=='A2'
                box+=(2*NNN-x_i[i])*r_ar[i]
            sum_C_R.append(round(box,5))
    else:
        for i in range(len(r_ar)):
            box=0
            #print(r_arr[i])
            for j in range(len(r_ar[i])):
                if par=='A1':
                    box+=x_i[j]*r_ar[i][j]
                    #print(' ',x_i[j],'*',r_arr2[i][j],end='')
                else: #par=='A2'
                    box+=(2*NNN-x_i[j])*r_ar[i][j]
                    #print(' ',2*NNN-x_i[j],'*',r_arr2[i][j],end='')
            sum_C_R.append(np.round(box,5))
    #print('\n',len(Mu),len(sum_C_R))
    print(par+':',mu,'\n',sum_C_R)
    print('\n\n')
    return sum_C_R
 
x_i=[x for x in range(0,11)]
#-------------------------------------------------------------

#fig, ax = plt.subplots(figsize=(12,9))
#"Scatter: dependence of A1 and A2 from Mu (∑(i=0..10)*r_j(i)"
import matplotlib.pyplot as plt
import pylab
import matplotlib.ticker
def point_graph(Mu,sum_C_R,sum_C_R2,sum_C_R3,sum_C_R4,f_size=(8,8),add=''):
    rcParams['figure.dpi']=150
    fig = plt.figure(figsize=f_size)
    main_target=catalogmaker(name='Graphics')
    
    ax1 = fig.add_subplot(211)
    ax1.set_title(u'Две области слиплись')
    ax2 = fig.add_subplot(212)

    ax1 = plt.subplot2grid((2, 1), (0, 0))
    ax1.set_title(u'Scatter: dependence of A1 from Mu (∑(i=0..10)*r_j(i)')
    ax2 = plt.subplot2grid((2, 1), (1, 0))
    ax2.set_title(u'Scatter: dependence of A2 from Mu (∑(i=0..10)*r_j(i)')
    
    ax1.scatter(x=Mu, y=sum_C_R, marker='o', c='r', edgecolor='b')
    ax1.scatter(x=Mu, y=sum_C_R3, marker='p', c='yellow', edgecolor='black')
    # Создаем экземпляр класса, который будет отвечать за расположение меток
    locator = matplotlib.ticker.MultipleLocator (base=0.1)#plt.xticks(np.arange(min(Mu), max(Mu)+0.1, 0.1))
    # Установим локатор для главных меток
    ax1.xaxis.set_major_locator (locator)

    ax2.scatter(x=Mu, y=sum_C_R2, marker='o', c='g', edgecolor='b')
    ax2.scatter(x=Mu, y=sum_C_R4, marker='p', c='purple', edgecolor='black')
    locator2 = matplotlib.ticker.MultipleLocator (base=0.1)
    ax2.xaxis.set_major_locator (locator2)

    ax1.set_xlabel('$Mu$')
    ax1.set_ylabel('$Avarage Number of A1 $')
    ax2.set_xlabel('$Mu$')
    ax2.set_ylabel('$Avarage Number of A2 $')

    for ax in fig.axes:
        ax.grid(True)

    plt.tight_layout(h_pad = 1)
    #save('pic_7_2_1', fmt='png')
    #save('pic_7_2_1', fmt='pdf')
    #plt.show()
    fig.savefig(main_target+'individs scatter_'+add)
    
#point_graph(Mu,sum_C_R,sum_C_R2,add='1_2')
#=============================================================

# # 3.+ Построить графики зависимости стационарных вероятностей при μ(1,2)≠μ(2,1) 
#инициализация матриц и поискстационарного распределения 
Mu= [0.1,0.3,0.5,0.8,1]
#def creator_mart_liambd(mu_12,mu_21=-1,NNN=5):
matrrr1=creator_mart_liambd(mu_12=Mu[0],mu_21=Mu[2])#;matr_printer(matrrr1,Mu[0],Mu[2]) 
matrrr2=creator_mart_liambd(mu_12=Mu[1],mu_21=Mu[2])#;matr_printer(matrrr2,Mu[1],Mu[2])
matrrr3=creator_mart_liambd(mu_12=Mu[2],mu_21=Mu[2])#;matr_printer(matrrr3,Mu[2],Mu[2])
matrrr4=creator_mart_liambd(mu_12=Mu[3],mu_21=Mu[2])#;matr_printer(matrrr4,Mu[3],Mu[2])
matrrr5=creator_mart_liambd(mu_12=Mu[4],mu_21=Mu[2])#;matr_printer(matrrr5,Mu[4],Mu[2])
'''=========================================='''        
#graphviz_plot(matrrr1,name_infile='fsm1',name_outfile='test_3_plus_gv',plottt=1)
'''=========================================='''
r_arr2=[[0.05072954],[0.0469718 ],[0.0142593 ],[0.00814817],[0.00679014],[0.00760496],
[0.0112346 ],[0.02246919],[0.06654338],[0.38024791],[0.38500101]]#0.8 
#присваиваем по началу произвольный массив из 11 элементов сумма которых 1, каждый <= 1, 
#который в последсвии будет состоять из 5 массив с 11 эл-в из прелдельных вероятносетей 

#R_stationar_finder(matrrr,mu_12,mu_21=-1,NNN=5):  
#matrrr(i)_transpose=np.transpose(matrrr(i)).copy()#Mu[(i-1)]
#r_arr2.append(R_stationar_finder(matrrr(i)_transpose,mu_12=Mu[(i-1)],mu_21=Mu[2]))
#Graph_Ri_i(r_arr2,Mu,Mu[2],add='not_simm')
#===============================================================

# # 6. Построить динамическую систему по 1 из матриц,которая удовлетворяет μ(1,2)≠μ(2,1) 
# ## 6.1 Ф-ция для вычисления значения усечения (численн. метод - дихотомия)  
# #произведем усечение, найдем границы интервала усечени А и В
# #scipy.stats.gamma()
# #np.random.gamma()#shape = k scale-theta
# #a_bord : x=[0, sscales * aas]
# #b_bord : x=[op=sscales * aas, sscales * aas * 2.5]
def F_calc(xxx,aa,sscale,Bord_a=0,Bord_b=0): 
    F=round(gamma.cdf(xxx,a=aa,scale=sscale),6)
    #print(F); x=[0,-> ~ 10] if scale =1/abs(liambda) 
    #a=k =2 or 3; scale = 1/abs(liambda) so it`s 1/|mart[i][i]|
    if Bord_a!=0:
        return round(F - Bord_a,6)
    elif Bord_b!=0:
        return round(F - Bord_b,6)
    else:
        return F

def dihotomy(x,aa=2,sscale=1/1,eps=0.01,bord='A',bordsc=0):#x=[a,b] a<b ,A_bord=0,B_bord=0
    x_1=[]
    for i in x:
        x_1.append(i)
    #tmp=[]
    counter=0
    while abs(x_1[1]-x_1[0])>eps:
            x_finalle=(x_1[0]+x_1[1])/2
            #print(counter,x_finalle,x_1)
            counter+=1
            if bord=='A':
                if F_calc(x_1[0],aa,sscale,Bord_a=bordsc)*F_calc(x_finalle,aa,sscale,Bord_a=bordsc)<0:
                    x_1[0]=x_1[0]
                    x_1[1]=x_finalle

                x_finalle=(x_1[0]+x_1[1])/2
                if F_calc(x_1[0],aa,sscale,Bord_a=bordsc)*F_calc(x_finalle,aa,sscale,Bord_a=bordsc)>0:
                    x_1[0]=x_finalle
                    x_1[1]=x_1[1]    
            else:#bord=='B':
                if F_calc(x_1[0],aa,sscale,Bord_b=bordsc)*F_calc(x_finalle,aa,sscale,Bord_b=bordsc)<0:
                    x_1[0]=x_1[0]
                    x_1[1]=x_finalle

                x_finalle=(x_1[0]+x_1[1])/2
                if F_calc(x_1[0],aa,sscale,Bord_b=bordsc)*F_calc(x_finalle,aa,sscale,Bord_b=bordsc)>0:
                    x_1[0]=x_finalle
                    x_1[1]=x_1[1]
    return x_finalle

# ## 6.2 Ф-ция поиска состояний для перехода и генерации времени 
def find_numb_OF_el(matr_l,i_state):#ищем на сколько и как делить отрезки,
    massiv=[]
    for i in range(len(matr_l)):
        if matr_l[i_state][i]>0:
            massiv.append([i,matr_l[i_state][i]])
    return massiv#ф-ция возвращает массив из !=0 сост, их индксы и лямбды
#--------------------------------------------------------------
def liambda_disturb_generator(matr,start_i,next_i,disturb_n='expon',kkk=1):#генерация времени прибывания с  усечением
    Alfa1_test=0.2#усечение
    Alfa2_test=1-Alfa1_test#для симметрии ограничений* 
    #liambda=abs(matr[start_i][start_i])#диагональный элемент
    liambda=abs(matr[start_i][next_i])
    if disturb_n=='expon':
        A_border=np.round(-np.log(1-Alfa1_test)/liambda,6)
        B_border=np.round(-np.log(1-Alfa2_test)/liambda,6)
        while True:
            time_generate=np.round(np.random.exponential(1/liambda),6)#depends of start_i : liambda=abs(matrix_to_r[start_i][start_i])
            if (A_border<=time_generate) and (B_border>=time_generate):
                break  
    else:#gamma 2k or 3k
        aas=kkk
        sscales=1/liambda
        op=sscales * aas * 2.5 #sscales * aas - mat ozhidanie
        A_border=dihotomy(x=[0,sscales * aas],aa=kkk,sscale=1/liambda,eps=0.01,bord='A',bordsc=0.2)
        B_border=dihotomy(x=[sscales * aas,op],aa=kkk,sscale=1/liambda,eps=0.01,bord='B',bordsc=0.8)
        while True:
            time_generate=np.asscalar(np.round(gamma.rvs(kkk,scale=1/liambda,size=1),6))
            if (A_border<=time_generate) and (B_border>=time_generate):
                break   
    return time_generate

# ## 6.3 Главная ф-ция - цикл, формирующий частоты и происводящий стох. процесс
def table4(matr,r_arr,i1,epsilon,dop=0,disturb_name='expon',k=1,main_target=os.getcwd() + os.sep):
    if not os.path.exists(main_target+str(i1)+disturb_name): 
        os.mkdir(main_target+str(i1)+disturb_name)# создание каталога 
        print('Каталог успешно создан', str(i1))
    target=main_target + str(i1) + disturb_name + os.sep
    
    #в первую очередь формируем таблицу
    with open (target+'table_last_%s_%s.txt'%(i1,dop),'w') as outfile:
        outfile.write('l'+'    '+'time_i'+'    '+'C(l)'+'    '+'Tau_i'+'    '+' Delta_n'+'\n')
    '''0ой блок - подготовка директории и файла-таблицы''' 
#--------------------------------------------------------------
    #R_i_n=np.zeros(len(r_arr))
    R_i_n2=np.zeros(len(r_arr))#numb of hits in each state
    v_i_n=np.zeros(len(r_arr))
    v_i_n2=np.zeros(len(r_arr))#frequency in each state
    tau_i=np.zeros(len(r_arr))
    tau_i2=np.zeros(len(r_arr))#time-frequency in each state
    
    TTT_i=np.zeros(len(r_arr))#all time in states: 
    #TTT_i_prev=np.zeros(len(r_arr))#for 7-8 task 
  
    time_i=np.array([[0]])#time
    start_i=i1;before_i=i1
    n=0 
    Delta_n=1
    '''1ый блок - подготовка переменны: массив частот,время пребываний в состояниях''' 
#--------------------------------------------------------------
    t_mimim=10
    while 1:#Delta_n>epsilon:
        
        if(n>0):
            v_i_n=v_i_n2.copy()
            tau_i=tau_i2.copy()
            #R_i_n=R_i_n2.copy()
            #TTT_i_prev=TTT_i.copy()
            before_i=start_i
            if t_mimim>time_generate:
                t_mimim=time_generate

        box = np.random.rand(1)#генерим ранд число куда пойдем
        segment = find_numb_OF_el(matr,start_i)#ищем на сколько и как делить отрезки
        '''2ой блок - начало запуска цикла + усл.остановки + после 1ой итер. запоминаем прошлые данные'''
#--------------------------------------------------------------
        if n>1:
            time_i=np.append(time_i,[time_i[-2]+time_i[-1]],axis = 0)
        #!!!time_i=np.append(time_i,[time_i[n-1]+time_generate],axis = 0)
        '''3ий блок - заполняем временную шкалу,время наступления события, те время когда сгенерим новое время'''
#--------------------------------------------------------------        
        n+=1
        summa2=0
        for i in range(len(segment)):
            summa2+=segment[i][1]

        summator=0;step=0
        for i in range(len(segment)):
            summator+=segment[i][1]/summa2#заменить len(segment) на сумму segment[i][1] i=1....n
            #print('summator',summator)
            if box <summator:
                if segment[i][0]>start_i:
                    step=1
                else:
                    step=-1
                start_i = segment[i][0]#в какое сотояние попали либо перешли в i-1 либо i+1
                break
        '''4ый блок - получаем в какое сост перейдем, генер. числа - 2ой блок'''
#-------------------------------------------------------------- 
        if step==1:# те мы шагнули вправо -> start_i=start_i+1
            time_generate=liambda_disturb_generator(matr,start_i-1,next_i=start_i,disturb_n=disturb_name,kkk=k)#    
        elif step==-1:# те мы шагнули влево -> start_i=start_i-1
            time_generate=liambda_disturb_generator(matr,start_i+1,next_i=start_i,disturb_n=disturb_name,kkk=k)#
        box=np.array([time_generate])#;print('\n',time_i)
        time_i=np.append(time_i,[box],axis = 0)
        '''5ый блок - add new полученное время пребывания в новом состоянии, обновляем частоты состояний'''
#--------------------------------------------------------------       
        TTT_i[before_i]=TTT_i[before_i]+time_generate#пребывание системы в опредлен состоянии(сколько время были) =abs(time_i[-2]-t_i[-1])    
        for i in range(len(tau_i2)):
            tau_i2[i]=np.round(TTT_i[i]/(time_i[-2]+time_i[-1]),6)
                
        if n>2:#критерий остановки по tau_i    
            Delta_n=round(max(abs(tau_i2-tau_i)),6)
            #Delta_n=round(max(abs(v_i_n2-v_i_n)),6)#TTT_i_prev-TTT_i_curr
            #print('delta:',Delta_n,'\n',TTT_i,'\n',TTT_i_prev)
            if Delta_n<=epsilon :# or n==2500:
                n-=1;print('break')
                with open (target+'table_last_%s_%s.txt'%(i1,dop),'a') as outfile:
                    outfile.write('n  '+'    '+str(round(time_i[-2][0],6))+'  to '+str(start_i)+'    '+str(round(abs(time_i[-1][0]),6))+'    '+str(Delta_n)+'\n')
                break
        #в случае если критерий по tau_i не выполнен мы переходим (добавлям в число попаданий +1) в новое состояние 
        #и отттуда снова генерим время        
        R_i_n2[start_i]+=1
        for i in range(len(tau_i2)):
            #tau_i2[i]=np.round(TTT_i[i]/(time_i[-2]+time_i[-1]),5)
            v_i_n2[i]=round(R_i_n2[i]/n,6)

        with open (target+'table_last_%s_%s.txt'%(i1,dop),'a') as outfile:
            outfile.write(str(n)+'    '+str(round(time_i[-2][0],6))+'    '+str(start_i)+'    '+str(round(abs(time_i[-1][0]),6))+'    '+str(Delta_n)+'\n')
        '''7ой блок - считаем дельту инкрементируем счетчик'''
#--------------------------------------------------------------
    v_i_n_test=np.array(v_i_n2)#(answ_mass[0])
    tau_i_test=np.array(tau_i2)#(answ_mass[1])
    with open (target+'table_last_%s_%s.txt'%(i1,dop),'a') as outfile:
        for ii in range(len(tau_i_test)+1):
            if ii == 0:
                outfile.write('\n\n' + 'i' + '\t'+'r_arr' + '\t\t' + 'v_i_n' + '\t\t' + 'tau_i' + '\n')
            else:
                i=ii-1
                #print(i,'\t',r_arr[i],'\t',v_i_n_test[i],'\t',tau_i_test[i])
                outfile.write(str(i)+'\t'+str(r_arr[i])+'\t'+str(v_i_n_test[i])+'\t\t'+str(tau_i_test[i])+'\n')
    an1=0;an2=0;an0=0
    for i in range(len(r_arr)):
        an0+=i*r_arr[i]
        an1+=i*v_i_n_test[i]
        an2+=i*tau_i_test[i]
    with open (target+'table_last_%s_%s.txt'%(i1,dop),'a') as outfile:
            outfile.write('\n\nAverage number of individs A1:\n' + 'r_arr:'+str(np.round(an0,5))+'\nv_i_n_test:'+str(np.round(an1,5))+'\ntau_i_test:'+str(np.round(an2,5))+'\n')
            outfile.write('\nR_i_n2:'+str(R_i_n2)+'\tsum(R_i_n2):'+str(sum(R_i_n2))+'\nv_i_n2:'+str(v_i_n2)+'\tsum(v_i_n2):'+str(sum(v_i_n2))+'\ntau_i2:'+str(tau_i2)+'\tsum(tau_i2):'+str(sum(tau_i2))+'\n')
    #print(time_i[-2]+time_i[-1])
    return [v_i_n_test,tau_i_test,R_i_n2,TTT_i,n,(time_i[-2]+time_i[-1])]
#===============================================================

# ## 6.4 Ф-ция - запускающия глав ф-цию N раз. и усредняет рез-ты. 
def N_times(matrix_to_r,r_arr,i1,eps,disturb_name='expon',kk=1,numbs=10,main_target=os.getcwd() + os.sep):#make 10 experiments for start place
    Global10_v_i_n=np.zeros(len(r_arr))#avarage count-frequency
    Global10_tau_i=np.zeros(len(r_arr))#avarage time-frequency
    Global10_R_i=np.zeros(len(r_arr))
    Global10_T_i=np.zeros(len(r_arr))
    N_generations=0
    all_Time=0
    Max_v_i_differ=np.zeros(len(r_arr))
    Max_tau_i_differ=np.zeros(len(r_arr))
    '''1ый блок - проводим эксперимент из i1-ого состояния numbs-раз и берем средние пакозатели'''
    for ten in range(numbs):
        anser=table4(matrix_to_r,r_arr,i1,eps,dop=ten,disturb_name=disturb_name,k=kk,main_target=main_target) 
        for i in range(len(r_arr)):
            Global10_v_i_n[i]+=np.round(anser[0][i]/numbs,5)
            Global10_tau_i[i]+=np.round(anser[1][i]/numbs,5)
            Global10_R_i[i]+=np.round(anser[2][i]/numbs)
            Global10_T_i[i]+=np.round(anser[3][i]/numbs,5)
        N_generations+=np.round(anser[4]/numbs)
        all_Time+=np.round(anser[5]/numbs,5)
    '''2ой блок - берем максимаьную разность от теор.стац. распред и v_i & tau_i + пишем все в файлы'''
    Max_v_i_differ[i1]=np.round(np.max(np.absolute(np.array(r_arr)-Global10_v_i_n)),5)
    Max_tau_i_differ[i1]=np.round(np.max(np.absolute(np.array(r_arr)-Global10_tau_i)),5)
    
    print('\n\nstart i=',i1,'\tN_times=',numbs,'\ni','\t','r_arr','\t\t','v_i_n','\t\ttau_i')
    with open (main_target+'table_finish_%s_%s.txt'%(i1,disturb_name),'w') as outfile:
        for ii in range(len(Global10_tau_i)+1):
            if ii == 0:
                outfile.write('N_times\n'+'i'+'\t'+'r_arr'+'\t\t'+'v_i_n'+'\t'+'tau_i'+'\n')
            else:
                i=ii-1
                print(i,'\t',r_arr[i],'\t',np.round(Global10_v_i_n[i],5),'\t',np.round(Global10_tau_i[i],5))
                outfile.write(str(i)+'\t'+str(r_arr[i])+'\t'+str(np.round(Global10_v_i_n[i],5))+'\t'+str(np.round(Global10_tau_i[i],5))+'\n')
    '''3ий блок - считаем кол-во особей типа А1 по теор.стац. распред и v_i & tau_i (усредненные)'''
    print('\nN_avar=',N_generations)
    print('\nT_avar=',all_Time)
    an1=0;an2=0;an0=0
    for i in range(len(r_arr)):
        an0+=i*r_arr[i]
        an1+=i*Global10_v_i_n[i]
        an2+=i*Global10_tau_i[i]
    print('Average number of individs A1:')
    print('r_arr:',np.round(an0,5))
    print('v_i_n_test:',np.round(an1,5))
    print('tau_i_test:',np.round(an2,5),'\n\n')
    with open (main_target+'table_finish_%s_%s.txt'%(i1,disturb_name),'a') as outfile:
            outfile.write('N_avar='+str(N_generations)+'\nT_avar='+str(all_Time)+'\n\nAverage number of individs A1:\n' + 'r_arr:'+str(np.round(an0,5))+'\nv_i_n_test:'+str(np.round(an1,5))+'\ntau_i_test:'+str(np.round(an2,5))+'\n')
            outfile.write('Max_v_i_differ='+str(Max_v_i_differ[i1])+'\nMax_tau_i_differ:'+str(Max_tau_i_differ[i1])+'\n')
            outfile.write('R_i_n2:='+str(Global10_R_i)+'\nv_i_n2:'+str(Global10_v_i_n)+'\ntau_i2:'+str(Global10_tau_i)+'\n')
    with open (main_target+'table_finals_%s.txt'%(disturb_name),'a') as outfile:
            outfile.write(str(i1)+'\t'+str(round(N_generations))+'\t'+str(np.round(all_Time,5))+'\t'+str(Max_v_i_differ[i1])+'\t\t'+str(Max_tau_i_differ[i1])+'\t\t\t'+str(np.round(an0,5))+'\t'+str(np.round(an1,5))+'\t\t'+str(np.round(an2,5))+'\n')    
    return [an0,an1,an2,Global10_v_i_n,Global10_tau_i,N_generations,all_Time,Global10_R_i,Global10_T_i]
    #{A1_by_r:an0,A1_by_v:an1,A1_by_tau:an2,v_i_n:Global10_v_i_n,tau_i:Global10_tau_i}
#===============================================================
def pd_saver(name,dicttt,targ=os.getcwd() + os.sep):
    tmp=dicttt
    temp_pd=pd.DataFrame(tmp)
    temp_pd.to_csv(targ+name)
    return temp_pd

def take_pd(name='r_arr_',Mu_name=0.8,column='r_i',add_place=0):
    if add_place==0:
        folder='pd_tables'
    else:
        folder=add_place + os.sep + 'pd_tables'
    targ=os.getcwd() + os.sep + folder + os.sep
    if name =='expon' or name =='gamma_2k' or name =='gamma_3k':
        y_answ0=[]
        y_answ1=[]
        y_answ2=[]
        _name=name
        y=pd.read_csv(targ+'y_answ_'+_name+'_'+str(Mu_name),sep=',',index_col=0)
        #print(y)
        for i in range(NNN*2+1):
            y_answ0.append(y['y_answ0'][i])
            y_answ1.append(y['y_answ1'][i])
            y_answ2.append(y['y_answ2'][i])
        return [y_answ0,y_answ1,y_answ2]
    
    else: # if name =='r_arr_'
        r=pd.read_csv(targ+name+str(Mu_name),sep=',',index_col=0)
        r_arr=[]
        for i in range(NNN*2+1):
                r_arr.append(round(r[column][i],8))
        return r_arr
#===============================================================
# 6.0 задание начальных условий (выбор матрицы); Mu= [0.1,0.3,0.5,0.8,1];
# В качестве примера возьмем A_μ(1,2)="0.8" ;μ(2,1)="0.5"
Mu_name=Mu[3]#Mu[3]
Delta_n=1;N_min=0#пока
matrix_to_r=matrrr4.copy()

r_arr=r_arr2[3].copy()#тк уже нашли стационарные ранее
r_temp=[]
for i in range(len(r_arr)):
    r_temp.append(np.round(r_arr[i][0],8))
    
pd_target=catalogmaker(name='pd_tables')
R_arr=pd_saver(name='r_arr_'+str(Mu_name),dicttt={'r_i':r_temp},targ=pd_target)
del r_temp;#print(R_arr)#;print(r_arr,'\n')

simmilar=[matrix_to_r[0][1]]
for i in range(1,len(matrix_to_r)):
    if i!=len(matrix_to_r)-1:
        simmilar.append(matrix_to_r[i][i-1])
        simmilar.append(matrix_to_r[i][i+1])
    else:
        simmilar.append(matrix_to_r[len(matrix_to_r)-1][len(matrix_to_r)-2])        
simmilar.sort()#print('\n',simmilar,len(simmilar),'\n',simmilar[0::6])
simmilar=pd_saver(name='simmilar_'+str(Mu_name),dicttt={'simmilar':simmilar},targ=pd_target)
#print(simmilar)
print(matrix_to_r)
print(R_arr,'\n')

#5
_name='gamma_2k' #'expon'#'gamma_2k'#'gamma_3k'
r=2        #1#2#3
amount=2
epsilon=0.00001

main_target=catalogmaker(name=_name)
with open (main_target+'table_finals_%s.txt'%(_name),'w') as outfile:
        outfile.write('i'+'\t'+'N'+'\t'+'T'+'\t\t'+'Max_v_i_differ'+'\t'+'Max_tau_i_differ'+'\t'+'r_arr'+'\t\t'+'v_i_n'+'\t\t'+'tau_i'+'\n')    
          
y_answ0,y_answ1,y_answ2=[],[],[]
v_i_pd,tau_i_pd=[],[]
i_begin,i_0_10=[],[]

T_toFIX,N_toFIX=[],[]
R_i_aver,T_i_aver=[0 for i in range(11)],[0 for i in range(11)]
v_i_aver,tau_i_aver=[0 for i in range(11)],[0 for i in range(11)]
import time
start_time = time.clock()
for i in range(11):    
    #a=N_times(matrix_to_r,r_arr,i,eps=epsilon,disturb_name=_name,kk=r,numbs=amount,main_target=main_target)
    y_answ0.append(round(a[0][0],5))
    y_answ1.append(round(a[1],5))
    y_answ2.append(round(a[2],5))
    for j in range(len(a[3])):
        v_i_pd.append(round(a[3][j],5))
        tau_i_pd.append(round(a[4][j],5))
        i_begin.append(i)
        i_0_10.append(j)
        
        v_i_aver[j]+=round(a[3][j]/11,5)
        tau_i_aver[j]+=round(a[4][j]/11,5)
        R_i_aver[j]+=a[7][j]/11 #тут мы складываем все получившиеся усрредненные R_i из стартовых i (11штук )
        T_i_aver[j]+=a[8][j]/11 #и делим их на количесво этих стартовых положений, то есть находим ОБЩЕЕ среднее
    N_toFIX.append(a[5])
    T_toFIX.append(np.asscalar(np.round(a[6],5)))
print ("{:g} s".format(time.clock() - start_time))
#y_answ;print(v_i_aver,sum(v_i_aver),'\n',tau_i_aver,sum(tau_i_aver))
#===============================================================
print(pd_target)
my_series=pd_saver(name='Freac_test_table_'+_name+'_'+str(Mu_name),dicttt={'v_i':v_i_pd,'tau_i':tau_i_pd,'i':i_0_10,'i_begin':i_begin},targ=pd_target)
y_answ=pd_saver(name='y_answ_'+_name+'_'+str(Mu_name),dicttt={'y_answ0':y_answ0,'y_answ1':y_answ1,'y_answ2':y_answ2},targ=pd_target)
del v_i_pd; del tau_i_pd; del i_0_10; del i_begin 
del y_answ0; del y_answ1; del y_answ2

T_toFIX=pd_saver(name='toFIX_T_'+_name+'_'+str(Mu_name),dicttt={'T_toFIX':T_toFIX},targ=pd_target)
N_toFIX=pd_saver(name='toFIX_N_'+_name+'_'+str(Mu_name),dicttt={'N_toFIX':N_toFIX},targ=pd_target)
R_i_aver=pd_saver(name='R_i_aver_'+_name+'_'+str(Mu_name),dicttt={'R_i_avar':np.round(R_i_aver,5)},targ=pd_target)
T_i_aver=pd_saver(name='T_i_aver_'+_name+'_'+str(Mu_name),dicttt={'T_i_avar':np.round(T_i_aver,5)},targ=pd_target)
v_i_aver=pd_saver(name='v_i_aver_'+_name+'_'+str(Mu_name),dicttt={'v_i_aver':np.round(v_i_aver,5)},targ=pd_target)
tau_i_aver=pd_saver(name='tau_i_aver_'+_name+'_'+str(Mu_name),dicttt={'tau_i_aver':np.round(tau_i_aver,5)},targ=pd_target)
print(np.round([1.234, 2.54]))
#===============================================================
#извлекаем данные
pre_name1='R_i_aver_';pre_name2='T_i_aver_'
pre_name3='toFIX_N_';pre_name4='toFIX_T_'
pre_name5='v_i_aver_';pre_name6='tau_i_aver_'
_name1='expon_';_name2='gamma_2k_';_name3='gamma_3k_'
R_2=take_pd(name=pre_name1+_name2,Mu_name=0.8,column='R_i_avar')#общее среднее от усредненных   
R_3=take_pd(name=pre_name1+_name3,Mu_name=0.8,column='R_i_avar')#кол-во вхождений в каждое сост
T_i_2=take_pd(name=pre_name2+_name2,Mu_name=0.8,column='T_i_avar')#общее среднее от усредненных
T_i_3=take_pd(name=pre_name2+_name3,Mu_name=0.8,column='T_i_avar')#время пребывания в каждом 
N_2=take_pd(name=pre_name3+_name2,Mu_name=0.8,column='N_toFIX')#среднее  
N_3=take_pd(name=pre_name3+_name3,Mu_name=0.8,column='N_toFIX')#кол-во поколений до фиксации из каждого 
T_2=take_pd(name=pre_name4+_name2,Mu_name=0.8,column='T_toFIX')#среднее  
T_3=take_pd(name=pre_name4+_name3,Mu_name=0.8,column='T_toFIX')#время до фиксации из каждого 
v_i_2=take_pd(name=pre_name5+_name2,Mu_name=0.8,column='v_i_aver')#общее среднее от усредненных
v_i_3=take_pd(name=pre_name5+_name3,Mu_name=0.8,column='v_i_aver')#частота по попаданиям
tau_i_2=take_pd(name=pre_name6+_name2,Mu_name=0.8,column='tau_i_aver')#общее среднее от усредненных
tau_i_3=take_pd(name=pre_name6+_name3,Mu_name=0.8,column='tau_i_aver')#частота по времени
ans=take_pd(name=pre_name4+_name2,Mu_name=0.8,column='T_toFIX')

folder='pd_tables'
folder=os.getcwd() + os.sep + folder + os.sep
Mu=0.8
_name='gamma_3k'
my_ser=pd.read_csv(folder+'Freac_test_table_'+_name+'_'+str(Mu),sep=',',index_col=0)
#print(my_ser)
list1=[]
list2=[]
for j in range(0,11):
    list1.append(round(my_ser['v_i'][j],5))
    list2.append(round(my_ser['tau_i'][ j],5))
print(0,list1,'\n',list2,'\n')

for i in range(1,11):
    list1=[]
    list2=[]
    print('\n\n',i)
    for j in range(11):
        print(round(my_ser['v_i'][i*11 + j ],5),',',end='')
        list1.append(round(my_ser['v_i'][i*11 + j],5))
    print()
    for j in range(11):
        print(round(my_ser['tau_i'][i*11 + j ],5),',',end='')
        list2.append(round(my_ser['tau_i'][i*11 + j],5))
    

name='gamma_3k'
y_answ0,y_answ1,y_answ2=take_pd(name,Mu_name=0.8)[0],take_pd(name,Mu_name=0.8)[1],take_pd(name,Mu_name=0.8)[2]
for j in range(0,11):
    print(np.round(y_answ1[j],5))
print()
for j in range(0,11):
    print(np.round(y_answ2[j],5))

#avarage_graph(y_answ0,y_answ1,y_answ2,name='gamma_2k')
for j in range(0,11):
    print(T_2[j])#N_2[j])T_2
print()
for j in range(0,11):
    print(T_3[j])#N_3[j])T_3

sum_C_R=average_funck(v_i_2,x_i,Mu,par='A1',len_ar=1)
sum_C_R=average_funck(v_i_3,x_i,Mu,par='A1',len_ar=1)
sum_C_R=average_funck(tau_i_2,x_i,Mu,par='A1',len_ar=1)
sum_C_R=average_funck(tau_i_3,x_i,Mu,par='A1',len_ar=1)
print(sum(v_i_2),'\n',sum(tau_i_2),'\n\n',sum(v_i_3),'\n',sum(tau_i_3),'\n')
print()
#===============================================================
r_erlang2=[]
r_erlang3=[]
names=['01','03','05','08','10']
Mu= [0.1,0.3,0.5,0.8,1]
type(names.index(names[0]))

for i in names:
    k2=take_pd(name='tau_i_aver_gamma_2k_',Mu_name=Mu[names.index(i)],column='tau_i_aver',add_place=i)
    r_erlang2.append(k2)
    #k3=take_pd(name='tau_i_aver_gamma_3k_',Mu_name=Mu[names.index(i)],column='tau_i_aver',add_place=i)
    #r_erlang3.append(k2)
    print(i,k2)
print('\n\n')
for i in names:
    k3=take_pd(name='tau_i_aver_gamma_3k_',Mu_name=Mu[names.index(i)],column='tau_i_aver',add_place=i)
    r_erlang3.append(k3)
    print(i,k3)

print('r_erlang2:')
sum_C_R_k2=average_funck(r_erlang2,x_i,Mu,par='A1')
sum_C_R2_k2=average_funck(r_erlang2,x_i,Mu,par='A2')
print('r_erlang3:')
sum_C_R_k3=average_funck(r_erlang3,x_i,Mu,par='A1')
sum_C_R2_k3=average_funck(r_erlang3,x_i,Mu,par='A2')

#Graph_Ri_i(r_erlang2,Mu,Mu[2],add='not_simm_k2')
#Graph_Ri_i(r_erlang3,Mu,Mu[2],add='not_simm_k3')
#point_graph(Mu,sum_C_R,sum_C_R2,add='2_2')
point_graph(Mu,sum_C_R_k2,sum_C_R2_k2,sum_C_R_k3,sum_C_R2_k3,add='4_2')