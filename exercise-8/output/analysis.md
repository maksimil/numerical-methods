# Надёжность

Данные для анализа в `part332.csv`.

**Одноэтапный метод (Метод Эйлера).** Хорошо оценивает погрешность.

**Двухэтапный Метод ($c_2=1/15$).** Недооценивает погрешность в <5 раз.

**Трёхэтапный Метод (Вторая Схема).** Недооценивает погрешность в <5 раз.

**Чётырёхэтапный Метод (Классический).** Переоценивает погрешность в >5 раз.

# Экономичность

Данные для анализа в `part333.csv`.

Для точностей до $10^{-2}$ **Метод Эйлера** требует примерно в 2 раза меньше
вызовов функции, чем остальные методы. Однако его число вызовов растёт быстрее
с ростом требуемой точности (он требует в 10 раз больше вызовов на точности
$10^{-7}$).

Остальные методы требуют примерно одного числа числа вызовов функции на
точностях до $10^{-5}$. **Двухэтапный Метод ($c_2=1/15$)** начинает требовать
в два раза больше вызовов на точностях выше.
