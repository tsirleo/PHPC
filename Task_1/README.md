# Практическое задание «Расписание сети сортировки» 

## Описание условия  

**Разработать последовательную программу вычисления:**

1. расписания сети сортировки, числа 
2. использованных компараторов и числа тактов, необходимых для её срабатывания при 
3. выполнении на n процессорах.  

Число тактов сортировки при параллельной обработке не должно превышать числа тактов, затрачиваемых четно-нечетной сортировкой Бетчера. 

Параметр командной строки запуска: n, где n>=1 – количество элементов в упорядочиваемом массиве, элементы которого расположены на строках с номерами [0…n-1]. 

***Формат команды запуска***: *bsort n* 

***Требуется***: 

1. вывести в файл стандартного вывода расписание и его характеристики в представленном далее формате; 
2. обеспечить возможность вычисления сети сортировки для числа элементов 1<=n<=10000; 
3. предусмотреть полную проверку правильности сети сортировки для значений числа сортируемых элементов 1<=n<=24; 
4. представить краткий отчет удовлетворяющий указанным далее требованиям. 

**Формат файла результата:** 

*Начало файла результата* 

*n 0 0* 

$cu_0$  $cd_0$ 

$cu_1$  $cd_1$ 

*…* 

$cu_{n_comp-1}$  $cd_{n_comp-1}$

*n_comp* 

*n_tact* 

*Конец файла результата* 

Где: 

*n 0 0*       – число сортируемых элементов, ноль, ноль. 

$cu_i$  $cd_i$   – номера строк, соединяемых i-м компаратором сравнения перестановки. 

*n_comp*  – число компараторов 

*n_tact*  – число тактов сети сортировки 


## Описание метода проверки 

Компиляция программы осуществляется командой: 
```shell
g++  bsort.cpp -o bsort
```

Запуск программы составления расписания для определенного числа n: 
```shell
./bsort n
```

Расписание выводится в заданном формате в файл “***schedule.txt***”. 

Запуск полной проверку правильности сети сортировки для количества сортируемых элементов ∈ [1, n]:
```shell
./bsort n -t
```

**Тестирование проводилось при помощи 0–1 принципа.**

0–1 принцип: если сеть сортирует все последовательности из нулей и единиц, то сеть является сортирующей. Необходимо перебрать все перестановки из n элементов, состоящие из нулей и единиц, пропустить их через сеть и проверить, что они корректно отсортированы. 

Была реализована функция ***generateBinaryArrays***, которая генерирует все перестановки из 0 и 1, длины n. Также была реализована функция ***testSchedule***, которая вызывается при тестировании и итеративно производит генерацию перестановок для каждого i ∈ [1, n], после чего составляется расписание и производится сортировку, а также вывод в файл.  

В результате работы тестовой программы создаются три файла: 

- “***input\_arrays.txt***” – вывод исходных массивов для каждого i ∈ [1, n]. 
- “***output\_arrays.txt***” – вывод отсортированных массивов для каждого i ∈ [1, n]. 
- “***comparators.txt***” – вывод расписания,  количества компараторов и  количества тактов для каждого i ∈ [1, n]. 
