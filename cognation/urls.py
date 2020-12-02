from django.conf.urls import url
from django.conf.urls.static import static
from django.conf import settings

from . import views

app_name = 'cognation'
urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^calculate/$', views.calculate, name='calculate'),
    url(r'^save_allele/', views.save_allele, name='save_allele')
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

