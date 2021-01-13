# paternity_genomed

<h1>1. Описание</h1>
<p>Данное приложение определяет вероятность родства в различных случаях. Полное описание формул, используемых в каждом случае, находится <a href='https://github.com/VieleSchritte/Genomed__software/tree/master/references'>в этих источниках</a>.</p>
<h2>Доступные случаи определения родства перечислены ниже.</h3>
<h3>1.1. Случаи с одним предполагаемым родителем:</h4>
<ul>
 <li>
  <h4>1.1.1. Отцовство/материнство для двух участников (отец и ребенок)</h4>
  <p>Исследованию подвергаются биологические образцы проверяемого ребенка и одного предполагаемого родителя. Биологический образец второго родителя для исследования недоступен</p>
 </li>
 <li>
  <h4>1.1.2. Отцовство/материнство для двух детей и одного предполагаемого родителя</h4>
  <p>Исследованию подвергают биологически образцы проверяемых детей и одного предполагаемого родителя. Биологическое родство между проверяемыми детьми достоверно известно. Биологический образец второго родителя для исследования недоступен</p>
 </li>
 <li>
  <h4>1.1.3. Отцовство/материнство для трех детей и одного предполагаемого родителя</h4>
  <p>Исследованию подвергают биологические образцы проверяемых детей и одного предполагаемого родителя. Биологическое родство между проверяемыми детьми достоверно известно. Биологический образец второго родителя для исследования недоступен</p>
 </li>
 <li>
  <h4>1.1.4. Отцовство/материнство для предполагаемого родителя, предполагаемого брата/сестры и проверяемого лица (родство между предполагаемыми родителем и братом/сестрой достоверно известно)</h4>
  <p>Исследованию подвергают биологические образцы проверяемого лица, одного предполагаемого родителя и одного предполагаемого брата или сестры. Предполагается, что проверяемое лицо и предполагаемый брат (сестра) имеют общих родителей, биологический образец одного из которых для исследования недоступен</p>
 </li>
</ul>
<h3>1.2. Случаи с предполагаемым и известным родителями:</h3>
<ul>
 <li>
  <h4>1.2.1. Один ребенок, известный и предполагаемый родители</h4>
  <p>Исследованию подвергают биологические образцы проверяемого ребенка, одного известного родителя и одного предполагаемого родителя. Биологическое родство между проверяемым ребенком и известным родителем достоверно известна</p>
 </li>
 <li>
  <h4>1.2.2. Двое детей, известный и предполагаемый родители</h4>
  <p>Исследованию подвергают биологические образцы проверяемых детей и одного предполагаемого родителя. Биологическое родство между проверяемыми детьми достоверно известно. Биологический образец второго родителя для исследования недоступен</p>
 </li>
 <li>
  <h4>1.2.3. Трое детей, известный и предполагаемый родители</h4>
  <p>Исследованию подвергают биологические образцы проверяемых детей, одного известного и одного предполагаемого родителя. Биологическое родство между проверяемыми детьми и известным родителем достоверно известно</p>
 </li>
</ul>
<h3>1.3. Случаи с родительской парой (достоверно известно, что это именно родительская пара, а не отдельные лица, с которыми проверяется родство):</h3>
<ul>
 <li>
  <h4>1.3.1. Один ребенок и родительская пара</h4>
  <p>Исследованию подвергают биологические образцы проверяемого ребенка и двух предполагаемых родителей. Достоверно известно, что предполагаемые родители являются родительской парой, а не отдельными лицами, с которыми проверяется родство</p>
 </li>
 <li>
  <h4>1.3.2. Двое детей и родительская пара</h4>
  <p>Исследованию повергают биологические образцы проверяемых детей и двух предполагаемых родителей. Достоверно известно, что предполагаемые родители являются родительской парой, а не отдельными лицами, с которыми проверяется родство. Также достоверно известно биологическое родство между проверяемыми детьми</p>
 </li>
 <li>
  <h4>1.3.3. Трое детей и родительская пара</h4>
  <p>Исследованию подвергают биологические образцы проверяемых детей и двух предполагаемых родителей. Достоверно известно, что предполагаемые родители являются родительской парой, а не отдельными лицами, с которыми проверяется родство. Также достоверно известно биологическое родство между проверяемыми детьми</p>
 </li>
</ul>
<h3>1.4. Случаи родства братьев и сестер:</h3>
<ul>
 <li>
  <h4>1.4.1. Родные братья/сестры</h4>
  <p>Исследованию подвергают биологические образцы проверяемого лица и одного предполагаемого брата или сестры. Предполагается, что исследуемые лица имеют общих родителей, биологические образцы которых для исследования недоступны</p>
 </li>
 <li>
  <h4>1.4.2. Проверяемое лицо и двое братьев и сестер</h4>
  <p>Исследованию подвергают биологические образцы проверяемого лица и двух предполагаемых братьев или сестер. Предполагается, что исследуемые лица имеют общих родителей, биологические образцы которых для исследования недоступны. Экспертная оценка более достоверна, чем при сравнении гаплотипов исследуемого лица и одного брата или сестры</p>
 </li>
 <li>
  <h4>1.4.3. Сводные/единокровные братья/сестры</h4>
  <p>Исследованию подвергают биологические образцы лиц, имеющих одного общего родителя. Это могут быть единокровные братья/сестры, при условии, что у них общий отец или единоутробные - если у них общая мать</p>
 </li>
 <li>
  <h4>1.4.4. Двоюродные братья/сестры</h4>
  <p>Исследованию подвергают биологические образцы проверяемого лица и предполагаемого двоюродного брата или сестры. Предполагается, что исследуемые лица имеют общих бабушку и дедушку, а один из родителей проверяемого лица является родным братом/сестрой одному из родителей второго участника. Биологические образцы родственников участников для исследования недоступны</p>
 </li>
</ul>
<h3>1.5. Случаи родства c участием бабушек и дедушек:</h3>
<ul>
 <li>
  <h4>1.5.1. Предполагаемая бабушка/дедушка и внук/внучка</h4>
  <p>Исследованию подвергают биологические образцы проверяемого лица и предполагаемой бабушки(дедушки). Предполагается, что один из родителей проверяемого лица является дочерью или сыном предполагаемой бабушки (дедушки). Биологические образцы предполагаемых родителей проверяемого лица для исследования недоступны</p>
 </li>
 <li>
  <h4>1.5.2. Бабушка, дедушка, внук/внучка</h4>
  <p>Исследованию подвергают биологические образцы проверяемого лица и предполагаемых бабушки и дедушки. Предполагается, что один из родителей проверяемого лица является дочерью или сыном предполагаемых бабушки и дедушки. Биологические образцы предполагаемых родителей проверяемого лица для исследования недоступны.</p>
 </li>
 <li>
  <h4>1.5.3. Проверяемое лицо, предполагаемый родитель, бабушка/дедушка (предполагаемый родитель является потомком бабушки/дедушки)</h4>
  <p>Исследованию подвергают биологические образцы проверяемого ребенка и его бабушки(дедушки) и одного предполагаемого родителя ребенка. Предполагается, что предполагаемый родитель ребенка является потомком бабушки (дедушки). Биологическое родство между проверяемым ребенком и бабушкой (дедушкой) достоверно известно. Биологический образец второго родителя ребенка для исследования недоступен</p>
 </li>
 <li>
  <h4>1.5.4. Проверяемое лицо, предполагаемый родитель, бабушка/дедушка (предполагаемый родитель НЕ является потомком бабушки/дедушки)</h4>
  <p>Исследованию подвергают биологические образцы проверяемого лица, предполагаемого родителя и предполагаемой бабушки/дедушки. Предполагается, что одним из родителей проверяемого лица является предполагаемый родитель, а другим - дочь или сын предполагаемой бабушки/дедушки. Биологический образец этого второго родителя для исследования недоступен</p>
 </li>
 <li>
  <h4>1.5.5. Проверяемое лицо, известный и предполагаемый родители, предполагаемая бабушка/дедушка</h4>
  <p>Исследованию подвергают биологические образцы проверяемого лица и предполагаемой бабушки(дедушки). Предполагается, что один из родителей проверяемого лица является дочерью или сыном предполагаемой бабушки (дедушки). Биологические образцы предполагаемых родителей проверяемого лица для исследования недоступны</p>
 </li>
 <li>
  <h4>1.5.6. Проверяемое лицо, предполагаемый родитель, предполагаемые бабушка и дедушка</h4>
  <p>Исследованию подвергают биологические образцы проверяемого лица, известного родителя, предполагаемой бабушки и предполагаемого дедушки. Предполагается, что одним из родителей проверяемого лица является известный родитель, а другим - дочь или сын предполагаемой бабушки и дедушки. Биологический образец этого второго родителя для исследования недоступен</p>
 </li>
</ul>
<h3>1.6. Прочие случаи:</h3>
<ul>
 <li>
  <h4>1.6.1. Проверяемое лицо и предполагаемый дядя/тетя</h4>
  <p>Исследованию подвергают биологические образцы проверяемого лица и дяди/тети. Предполагается, что один из родителей проверяемого лица является родным братом или сестрой второго учасника. Биологический образец этого родителя для исследования недоступен</p>
 </li>
 <li>
  <h4>1.6.2. Проверяемое лицо и предполагаемая бабушка/дедушка (обсчет с помощью IBD-индексов, подробнее - <a href='https://github.com/VieleSchritte/Genomed__software/tree/master/references'>в источниках</a>)</h4>
  <p>Исследованию подвергают биологические образцы проверяемого лица и предполагаемой бабушки(дедушки). Предполагается, что один из родителей проверяемого лица является дочерью или сыном предполагаемой бабушки (дедушки). Биологические образцы предполагаемых родителей проверяемого лица для исследования недоступны</p>
 </li>
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
  <p>Приложение поддерживает возможность ввода дробных значений аллелей. В качестве разделителя принимаются как запятые, так и точки. При вводе запятой в качестве разделяющего символа приложение проведет расчет, заменив в таблице результатов запятую на точку:</p>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/comma_delimiting.png">
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/comma_del_res.png">
 </li>
 <li>В данном случае всего два участника, но в случаях, предполагающих большее количество лиц, нужно будет нажать кнопку "Далее" для дальнейшего ввода гаплотипов. Если гаплотипы введены неправильно, но кнопка "Далее" уже нажата, можно вернуться на предыдущий ввод по кнопке "Назад". При этом все гаплотипы сохранятся. Также гаплотипы сохраняются, если по кнопкам "Назад" был совершен переход к описанию. Гаплотипы удаляются, если с описания совершен переход назад к выбору родства, т.к. при выборе другого вида родства необходимо корректно передать данные на сервер.</li>
 <li>
  <p>По кнопке "Рассчитать" передаем данные на сервер для обсчета и получаем отчет:</p>
  <img src="https://github.com/VieleSchritte/Genomed__software/blob/master/readme_files/final_table.png">
 </li>
 <li>В отчете представлены гаплотипы участников по каждому локусу и вероятность родства, обозначенная аббревиатурой LR (Likehood Ratio). Также в конце таблицы представлены значения CPI (Combined Paternity Index - комбинированного индекса отцовства, полученного произведением всех значений LR в таблице) и вероятности отцовства, вычисляемой по формуле: P = CPI / (1 + CPI)</li>
 <li>
  <p>Для скачивания таблицы в формате .xlsx необходимо нажать "Сохранить результаты". По кнопке "К выбору родства" можно при необходимости перейти на стартовую страницу.</p>
  <p>Если получена таблица с результатами расчета, то все сделано правильно.</p>
 </li>
</ol>
<h1>3. Дополнение в случае гомозиготного состояния</h1>
<p>В случае гомозиготности по какому-либо локусу, кроме полоспецифичных (см. ниже), приложение поддерживает дополнение аллелей: при введении одного числа в финальной таблице оно автоматически проставляет через слэш то же число.</p>
<img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/homozygosity_enter.png">
<img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/homozygosity_output.png">
<h1>4. Ввод аллелей при работе с полоспецифичными локусами.</h1>
<p>В наборах, поддерживаемых данным приложением, используются следующие полоспецифичные локусы:</p>
<ul>
 <li>AMEL</li>
 <li>Yindel</li>
 <li>DYS391</li>
 <li>SRY</li>
</ul>
<p>Данные локусы выводятся в конце таблицы результатов и не участвуют в расчете. Все они, кроме амелогенина (AMEL), представлены одной копией, т.к. представлены только на Y-хромосоме, поэтому для них корректно выводить в финальную таблицу только один аллель. Таким образом, приложение проверяет все локусы, введенные для обоих участников, и, если находит среди них полоспецифичные, не добавляет через слэш еще один аллель, как в случае с гомозиготным состоянием:</p>
<img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/correct_enter_gender_specific.png">
<img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/correct_gender_specific_output.png">
<p>В случае амелогенина (AMEL) поддерживается дополнение гомозиготного состояния:</p>
<img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/AMEL_homozygous.png">
<img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/AMEL_homozygous_output.png">
<p></p>
<p>Случаи, в которых введено два аллеля, являются для полоспецифичных локусов ошибочными. Эти случаи обрабатываются в исключении, предусмотренном для слишком большого числа аллелей, которое обрабатывается через рендер страницы исключений:</p>
<img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/2_alleles_gender_spec.png">
<img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/2_alleles_gender_spec_exception.png">
<p>Подробнее об ошибках, требующих обработки исключений, - в следующем разделе.</p>
<h1>5. Обработка исключений</h1>
<p>Реализована обработка возможных ошибок ввода гаплотипов. Так, приложение подскажет, что было введено неправильно, и в каком локусе была допущена ошибка. Данные ошибки обрабатываются через рендер страницы исключений. Примеры возможных ошибок представлены ниже.</p>
<ol>
 <li>
  <h3>Допущена опечатка в названии локуса:</h3>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/wrong_locus.png">
  <p>В этом случае будет отображено, в каком локусе была допущена ошибка - выведется неправильно введенное название локуса.</p>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/wrong_locus_exception.png">
 </li>
 <li>
  <h3>При вводе длины допущена ошибка: символ, не являющийся числом или разделителем:</h3>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/wrong_symbol_allele.png">
  <p>При этом также будет выведен символ, который является ошибочным, локус и пара аллелей, в которых допущена ошибка:</p>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/wrong_symbol_allele_exception.png">
 </li>
 <li>
  <h3>Нет символов после разделителя или разделяющий символ стоит первым:</h3>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/del_first.png">
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/del_last.png">
  <p>На странице исключений происходит визуализация ошибки с указанием локуса:</p>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/del_first_exception.png">
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/del_last_exception.png">
 </li>
 <li>
  <h3>Аллели не введены</h3>
  <p>Ошибка покрывает случай, когда введен локус, но не введены соответствующие ему аллели. При этом могут быть введены пробельные символы:</p>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/empty_alleles.png">
  <p>В данном случае на странице исключений будет выведен локус, в котором это произошло: </p>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/empty_alleles_exception.png">
 </li>
 <li>
  <h3>Введено слишком большое число аллелей</h3>
  <p>Данный случай частично рассматривался в разделе 4 "Ввод аллелей при работе с полоспецифичными локусами". Для локусов, не являющихся полоспецифичными, неправильным является число аллелей, большее 2:</p>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/too_much_alleles_enter.png">
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/too_much_alleles_output.png">
 </li>
 <li>
  <h3>Неправильный формат строки</h3>
  <p>Если формат строки не соответствует введению локуса и аллелей, то приложение уведомит об этом:</p>
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/bullshit_enter.png">
  <img src="https://github.com/VieleSchritte/paternity_genomed/blob/master/readme_files/bullshit_exception.png">
 </li>
</ol>

<h1>6. Специфичный случай исключения</h1>
<p>Приложение поддерживает добавление нового значения аллеля. Разберем такое добавление на примере локуса TPOX:</p>
  <img src="https://github.com/VieleSchritte/Genomed__software/blob/master/readme_files/wrong_length_exception.png">
  <p>В этом случае в строке таблицы результатов, соответствующей этому локусу, будет выведено следующее:</p>
  <img src="https://github.com/VieleSchritte/Genomed__software/blob/master/readme_files/wrong_length_raw.png">
  <p>При нажатии на кнопку "ADD" появится окно для ручного ввода частоты:</p>
  <img src="https://github.com/VieleSchritte/Genomed__software/blob/master/readme_files/frequency_addition.png">
  <p>После нажатия на кнопку "Save" частота будет занесена в базу, а расчет - переделан с учетом введенной величины.</p>
