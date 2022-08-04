from django.shortcuts import render

def home(request):
    return render(request, 'djangopage/home.html')

# blog->templates ->blog ->template.html  (Logic about creating a template and blog folder(django convention))
