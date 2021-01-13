# paternity_genomed

<h1>1. Описание</h1>
<p>Данное приложение определяет вероятность родства в различных случаях. Полное описание формул, используемых в каждом случае, находится <a href='https://github.com/VieleSchritte/Genomed__software/tree/master/references'>в этих источниках</a>.</p>
<h3>1.1. Доступные случаи определения родства перечислены ниже.</h3>
<h4>1.1.1. Случаи с одним предполагаемым родителем:</h4>
<ul>
 <li>1.1.1.1. Отцовство/материнство для двух участников (отец и ребенок)</li>
 <li>1.1.1.2. Отцовство/материнство для двух детей и одного предполагаемого родителя</li>
 <li>1.1.1.3. Отцовство/материнство для трех детей и одного предполагаемого родителя</li>
 <li>1.1.1.4. Отцовство/материнство для предполагаемого родителя, предполагаемого брата/сестры и проверяемого лица (родство между предполагаемыми родителем и братом/сестрой достоверно известно)</li>
</ul>
<h4>1.1.2. Случаи с предполагаемым и известным родителями:</h4>
<ul>
 <li>1.1.2.1. Один ребенок, известный и предполагаемый родители</li>
 <li>1.1.2.2. Двое детей, известный и предполагаемый родители</li>
 <li>1.1.2.3. Трое детей, известный и предполагаемый родители</li>
</ul>
<h4>1.1.3. Случаи с родительской парой (достоверно известно, что это именно родительская пара, а не отдельные лица, с которыми проверяется родство):</h4>
<ul>
 <li>1.1.3.1. Один ребенок и родительская пара</li>
 <li>1.1.3.2. Двое детей и родительская пара</li>
 <li>1.1.3.3. Трое детей и родительская пара</li>
</ul>
<h4>1.1.4. Случаи родства братьев и сестер:</h4>
<ul>
 <li>1.1.4.1. Родные братья/сестры</li>
 <li>1.1.4.2. Проверяемое лицо и двое братьев и сестер</li>
 <li>1.1.4.3. Сводные/единокровные братья/сестры</li>
 <li>1.1.4.4. Двоюродные братья/сестры</li>
</ul>
<h4>1.1.5. Случаи родства c участием бабушек и дедушек:</h4>
<ul>
 <li>1.1.5.1. Предполагаемая бабушка/дедушка и внук/внучка</li>
 <li>1.1.5.2. Бабушка, дедушка, внук/внучка</li>
 <li>1.1.5.3. Проверяемое лицо, предполагаемый родитель, бабушка/дедушка (предполагаемый родитель является потомком бабушки/дедушки)</li>
 <li>1.1.5.4. Проверяемое лицо, предполагаемый родитель, бабушка/дедушка (предполагаемый родитель НЕ является потомком бабушки/дедушки)</li>
 <li>1.1.5.5. Проверяемое лицо, известный и предполагаемый родители, предполагаемая бабушка/дедушка (один из родителей является потомком бабушки/дедушки)</li>
 <li>1.1.5.6. Проверяемое лицо, предполагаемый родитель, предполагаемые бабушка и дедушка</li>
</ul>
<h4>1.1.6. Прочие случаи:</h4>
<ul>
 <li>1.1.6.1. Проверяемое лицо и предполагаемый дядя/тетя</li>
 <li>1.1.6.2. Проверяемое лицо и предполагаемая бабушка/дедушка (обсчет с помощью IBD-индексов, подробнее - <a href='https://github.com/VieleSchritte/Genomed__software/tree/master/references'>в источниках</a>)</li>
</ul>
<h1>2. Использование</h1>
<p>Рассмотрим использование данного сервиса на примере получения вероятности отцовства в случае, когда доступны гаплотипы проверяемого лица (ребенка) и предполагаемого родителя. Гаплотипы участников нужно будет взять <a href="https://github.com/VieleSchritte/Genomed__software/blob/master/references/trial_data_parent">отсюда.</a></p>
<ol>
 <li>
  <p>Заходим на сайт приложения: https://genomed-paternity.ru/ (логин и пароль находятся локально на компьютерах в лаборатории), в интерфейсе выбираем нужный тип родства:</p>
  <img src="https://github.com/VieleSchritte/Genomed__software/blob/master/readme_files/start_screen.png">
 </li>
 <li>
  <p>Для лучшего понимания того, с какими данными работает каждая формула, в каждый случай родства добавлено описание со схемой. Настоятельно рекомендуется ознакомиться с этим материалом, чтобы получить корректный результат.</p>
  <img src="https://github.com/VieleSchritte/Genomed__software/blob/master/readme_files/description_screen.png">
 </li>
 <li>По кнопке "ОК" переходим ко вводу гаплотипов</li>
 <li>
  <p>Для каждого гаплотипа предусмотрено отдельное поле ввода с описанием, так что перепутать гаплотипы не получится. Вводим каждый гаплотип в свое окно:</p>
  <img src="https://github.com/VieleSchritte/Genomed__software/blob/master/readme_files/genotypes_enter.png">
 </li>
 <li>В данном случае всего два участника, но в случаях, предполагающих большее количество лиц, нужно будет нажать кнопку "Далее" для дальнейшего ввода гаплотипов. Если гаплотипы введены неправильно, но кнопка "Далее" уже нажата, можно вернуться на предыдущий ввод по кнопке "Назад". При этом все гаплотипы сохранятся. Также гаплотипы сохраняются, если по кнопкам "Назад" был совершен переход к описанию. Гаплотипы удаляются, если с описания совершен переход назад к выбору родства, т.к. при выборе другого вида родства необходимо корректно передать данные на сервер.</li>
 <li>
  <p>По кнопке "Рассчитать" передаем данные на сервер для обсчета и получаем отчет:</p>
  <img src="https://github.com/VieleSchritte/Genomed__software/blob/master/readme_files/final_table.png">
 </li>
 <li>В отчете представлены гаплотипы участников по каждому локусу и вероятность родства, обозначенная аббревиатурой LR (Likehood Ratio). Также в конце таблицы представлены значения CPI (Combined Paternity Index - комбинированного индекса отцовства, полученного произведением всех значений LR в таблице) и вероятности отцовства, вычисляемой по формуле: P = CPI / (1 + CPI)</li>
 <li>Для скачивания таблицы в формате .xlsx необходимо нажать "Сохранить результаты". По кнопке "К выбору родства" можно при необходимости перейти на стартовую страницу</li>
</ol>
Если получена таблица с результатами расчета, то все сделано правильно.
<h1>3. Обработка исключений</h1>
<p>Реализована обработка возможных ошибок ввода гаплотипов. Так, приложение подскажет, что было введено неправильно, и в каком локусе была допущена ошибка. Примеры возможных ошибок представлены ниже.</p>
<ol>
 <li>
  <h3>Введен аллель, которого нет в базе данных - в данном случае - в локусе TPOX:</h3>
  <img src="https://github.com/VieleSchritte/Genomed__software/blob/master/readme_files/wrong_length_exception.png">
  <p>В этом случае в строке таблицы результатов, соответствующей этому локусу, будет выведено следующее:</p>
  <img src="https://github.com/VieleSchritte/Genomed__software/blob/master/readme_files/wrong_length_raw.png">
  <p>При нажатии на кнопку "ADD" появится окно для ручного ввода частоты:</p>
  <img src="https://github.com/VieleSchritte/Genomed__software/blob/master/readme_files/frequency_addition.png">
  <p>После нажатия на кнопку "Save" частота будет занесена в базу, а расчет - переделан с учетом введенной величины.</p>
 </li>
 <li>
  <h3>Допущена опечатка в названии локуса:</h3>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/wrong_locus.png">
  <p>В этом случае при отправке данных на сервер исключение будет обработано через рендер страницы исключений, на которой будет отображено, в каком локусе была допущена ошибка - выведется неправильно введенное название локуса.</p>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/wrong_locus_exception.png">
 </li>
 <li>
  <h3>При вводе длины допущена ошибка: символ, не являющийся числом, и не являющийся разделителем:</h3>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/wrong_symbol_allele.png">
  <p>Такие ошибки также обрабатываются через рендер страницы исключений. При этом будет выведен символ, который является ошибочным, локус и пара аллелей, в которых допущена ошибка:</p>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/wrong_symbol_allele_exception.png">
  <p>При вводе запятой в качестве разделяющего символа программа все равно проведет расчет, заменив в таблице результатов запятую на точку:</p>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/comma_delimiting.png">
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/comma_del_res.png">
 </li>
 
</ol>
