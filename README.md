  Данный контейнерный класс представляет собой структуру для математической работы с матрицами.
  В классе реализованы два конструктора. В первом один параметр по умолчанию (им заполняются все ячейки матрицы) и два беззнаковых числа – размеры матрицы, 
в данном конструкторе создается двумерный динамический массив. Второй – конструктор копирования. Прописан и деструктор, который очищает созданный 
двумерный динамический массив.

  Были перегружены следующие операторы:
    •	«*», который математически умножает матрицу на число или на другую матрицу.
    •	«/» – обратное действие к умножению матрицы на число.
    •	«+», который математически складывает матрицы.
    •	«-» – обратное действие к сумме матриц.
    •	«=» и другие операторы присваивания, совмещенные с вышеперечисленными операциями.
  Методы getRows() и getCols() возвращают соответственно количество строк и столбцов матрицы.
  Метод transposition() возвращает транспонированную матрицу.
  Было реализовано два метода расчета определителя, каждый из которых имеет свои преимущества и недостатки.
    •	Расчет определителя через алгебраическое дополнение. Минусом этого метода является сложность по времени, которая растет по экспоненте, 
    если размер квадратной разницы больше 2. Плюсом является его точность, поскольку потерь при делении чисел не возникает.
    •	Расчет определителя методом Гаусса. По результатам тестов, данный метод работает намного быстрее предыдущего, однако, точность 
    начинает падать при размерах матрицы выше 10. Его сложность равняется О(n3).
  Был реализован общий метод вычисления определителя, который, в зависимости от размерности матрицы вызывает один из двух вышеперечисленных 
методов (что не запрещает пользоваться публичными методами расчета определителя). Данный метод в первую очередь должен быть точным, и лишь во вторую быстрым. 
Поэтому для матрицы размерности от 3 до 9 включительно вызывается метод Гаусса, а во всех остальных случаях метод алгебраического дополнения.
  Самый главный метод – расчет обратной матрицы. Для него было реализовано три вспомогательных метода – расчета минора и алгебраического дополнения, 
а также транспонирование матрицы. Данный метод возвращает матрицу типа double (во избежание ошибок, поскольку элементы обратной матрицы, как правило, дробные числа).
  
  Для проведения тестов было реализовано 2 класса: TestMatrix и Timer. В первом классе написаны Unit тесты, второй создан для замера времени работы 
выбранного куска кода.
  Метод для сравнения точности алгоритма Гаусса принимает в качестве параметра размер матрицы и количество проверок. Логика работы данного метода следующая: 
создается матрица определенного размера, рассчитывается определитель через алгебраическое дополнение и метод Гаусса, потом идет сравнение двух полученных величин, 
допускаемая погрешность: ε = 10-7. Если разница двух величин больше чем ε, то данный метод возвращает false. И так определенное количество итераций. 
Если цикл завершается не раньше положенного, то данный метод возвращает true.
  Метод для сравнения времени работы принимает в качестве параметра размер матрицы и количество проверок, а возвращает массив из двух элементов – 
время расчета определителя методом алгебраического дополнения и методом Гаусса соответственно определенное количество раз.

  Пример использования класса представлен в конце файла main.cpp. Там создается матрица и выводится обратная к ней.
  
  Оператор << был перегружен исключительно для удобства тестирования, пользоваться им в для реального вывода матрицы не рекомендуется.
