function f=benchmark_func(x,func_num)
format long e;
global initial_flag;
   initial_flag = 0;
   [f,g,h]=CEC2017(x,func_num);
end
