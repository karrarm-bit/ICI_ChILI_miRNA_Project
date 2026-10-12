# دليل رفع المشروع والحصول على DOI

## الرفع

1. فك ضغط الملف. افتح مجلد `ICI_ChILI_miRNA_Project`.
2. أنشئ مستودعًا على GitHub، والاسم المقترح `ICI_ChILI_miRNA_Project`.
3. من `Add file > Upload files` ارفع **محتويات المجلد** مع الحفاظ على
   المجلدات الفرعية. لا ترفع ZIP وحده، ولا تجعل README داخل مجلد إضافي.
4. افحص الصفحة: يجب أن يظهر README وscripts وdata وresults وdocs في الجذر.
5. أكمل أسماء الطلبة والمشرف والفرع في المسودة. أكمل ملف
   `CITATION.cff.example` ثم غيّر اسمه إلى `CITATION.cff`، مع حذف كلمة
   REPLACE واستبدالها بالبيانات الحقيقية. لا تنشر بيانات مؤلفين غير مؤكدة.
6. اختر رخصة للكود بصفتكم أصحاب الحقوق. القالب الحالي لا يفرض رخصة.
7. شغّل `VERIFY_REPOSITORY.R` ثم وضع `review`، ويفضل `full` قبل الإصدار.
   هذه الحزمة فُحصت نحويًا وبفحص البيانات، لكن لم يُنفذ edgeR هنا.

GitHub يسمح برفع حتى 100 ملف في الدفعة من المتصفح وبحد 25 MiB للملف.
الحزمة معدّة ضمن هذين الحدين. هذه حدود الرفع بالمتصفح، وليست حدًا
لإجمالي الملفات التي يمكن أن يحتويها مستودع GitHub.

## DOI عبر Zenodo

1. سجّل الدخول إلى https://zenodo.org وربط حساب GitHub بحساب Zenodo.
2. افتح قائمة الحساب في Zenodo ثم GitHub، واضغط Sync now إذا لزم.
3. فعّل المستودع المقصود من المفتاح المقابل لاسمه.
4. بعد تفعيل الربط، أنشئ Release جديدًا في GitHub، مثل `v1.0.0`.
5. انتظر المعالجة ثم افتح سجل Zenodo الناتج وخذ DOI الحقيقي.
6. راجع أسماء المؤلفين والعنوان والرخصة والملفات في سجل Zenodo.
7. أضف DOI الحقيقي إلى README وCITATION.cff وإلى بحث التخرج.

لا تنشئ رقم DOI افتراضيًا. تجهيز ZIP لا يولد DOI؛ الإصدار والأرشفة هما
الخطوة المطلوبة. إذا سبق وجود DOI لهذا **المشروع نفسه**، استخدم الإصدار
التالي من السجل القائم بدل إنشاء مشروع مكرر. DOI مشروع آخر لا يُعاد
استخدامه لهذا البحث.

## كيف تتأكد أن الرفع مكتمل؟

بعد الرفع نزّل نسخة المستودع من `Code > Download ZIP`. فكها وقارن
الملفات مع `FILE_MANIFEST.csv`. يمكن التحقق من SHA256SUMS.txt باستخدام
`sha256sum -c SHA256SUMS.txt` على Linux/Git Bash، أو Get-FileHash على Windows.
ملفا manifest وchecksums لا يسجلان بصمتهما الذاتية لتجنب الاعتماد الدائري؛
يسجل manifest اسم SHA256SUMS.txt دون بصمة. نجاح تجهيز الحزمة هنا لا يثبت
اكتمال الرفع الذي ستجريه لاحقًا؛ المقارنة بعد التنزيل تثبت سلامة النسخة.

## ملاحظة علمية للمقيّم

النتيجة الرئيسية للجزء المتعلق بالـmiRNA هي مراجعة hairpin-associated:
65 feature، و12 عينة، وmiR-26b-5p فقط وفق FDR ومعيار fold change.
الأربع مرشحات السابقة تخص التحليل المجمّع المساند. تقاطع 17 جينًا لم
يظهر إثراءً إحصائيًا، ولا يثبت آلية تنظيم مباشرة بين مجموعتين مستقلتين.

المراجع الرسمية:
- https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository
- https://help.zenodo.org/docs/profile/linking-accounts/
- https://help.zenodo.org/docs/github/enable-repository/
- https://help.zenodo.org/docs/github/archive-software/github-upload/
- https://help.zenodo.org/docs/github/describe-software/
