# VKR
Moran model analysis

Тема ВКР:
"Расчёт характеристик генетической эволюции в стохастических моделях Морана".
=====================

Уточнение по модели Морана с мутацией:
-----------------------------------
* 1 диплоидная модель, N=5 (10 гаплоидных особей c аллелями A1 и A2);
* 2 состояние популяции в любой момент времени характеризуется
числом i аллелей A1 (i=0,...,10),
i=0 означает, что все особи имеют аллели A2,
i=10 означает, что все особи имеют аллели A1.
* 3 матрица интенсивностей переходов между состояниями:  
λ(0,1)=2Nμ(2,1), λ(0,0)=-2Nμ(2,1);  
λ(N,N-1)=2Nμ(1,2), λ(N,N)=-2Nμ(1,2);  
λ(i,i-1)=[i(2N-i)/2N]+μ(2,1)(2N-i), λ(i,i+1)=[i(2N-i)/2N]+μ(1,2)i;  
λ(i,i)=-μ(i,i-1)-μ(i,i+1) при 0<i<2N, где коэффициент λ(1,2) является интенсивностью мутации одной особи с аллелью A1 в особь с аллелью A2,
коэффициент λ(2,1) является интенсивностью мутации одной особи с аллелью A2 в особь с аллелью A1.

Задание:
* 1 нарисовать графы системы для поглощающей модели и без фиксирующих состояний;
* 2 составить систему уравнений Колмогорова для вероятностей состояний;
* 3 при значениях коэффициентов мутации μ=μ(1,2)=μ(2,1)=0,1 ; 0,3 ; 0,5 ; 0,8 ; 1
найти стационарные распределения и построить графики зависимости стационарных вероятностей r(i) от i=0,...,10;
* 4 построить график зависимости среднего числа
особей c аллелью A1 от ?=?(1,2)=?(2,1) (μ меняется от 0,1 до 1 с шагом 0,05);
* 5 построить график зависимости среднего числа особей c аллелью A2 от μ=μ(1,2)=μ(2,1)
* 6 произвести расчет характеристик модели для распределение Эрланга(к=2,3), показательного
* 7 построить графики соотношений характеристик
* 8 провести анализ полученных результатов, сделать выводы.
