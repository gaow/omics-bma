GSL_VERSION := 1.16
MOSEK_VERSION := 7.1.0.44
MOSEK_PATH := $(HOME)/.mosek
MOSEK_LIC := $(HOME)/.mosek.lic
swig_opts := -c++ -python -O -shadow -keyword -w-511 -w-509

install:
	python setup.py install

install_libs: igsl imosek

ieqtlbma:
	swig $(swig_opts) -o src/pyeqtlbma/pyeqtlbma_wrap.cxx src/pyeqtlbma/pyeqtlbma.i && \
	mv src/pyeqtlbma/pyeqtlbma.py src/pyeqtlbma/__init__.py

igsl:
	cd external/gsl && \
	./configure CFLAGS="-O3 -fPIC" --prefix=`pwd` && make && make install

imosek: $(MOSEK_LIC)
	cp -a external/mosek $(MOSEK_PATH)
	cd $(MOSEK_PATH)/7/tools/platform/linux64x86/python/3 && \
	python setup.py install
	@echo "\033[0;32m~~~~~~~~~~~~~~ IMPORTANT ~~~~~~~~~~~~~~~~"
	@echo 'Please make sure MOSEK environment is set'
	@echo 'For bash in Linux, run:'
	@echo 'echo "export PATH=$(MOSEK_PATH)/7/tools/platform/linux64x86/bin:\$PATH" >> ~/.bashrc'
	@echo 'echo "export LD_LIBRARY_PATH=$(MOSEK_PATH)/7/tools/platform/linux64x86/bin:\$LD_LIBRARY_PATH" >> ~/.bashrc'
	@echo 'echo "export MOSEKLM_LICENSE_FILE=$(MOSEK_LIC)" >> ~/.bashrc'
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\033[0m"

download: eqtlbma gsl mosek

eqtlbma:
	cd external && \
	git clone -b omics-bma https://github.com/gaow/eqtlbma.git

gsl:
	cd external && \
	rm -rf gsl* && \
	wget ftp://ftp.gnu.org/gnu/gsl/gsl-$(GSL_VERSION).tar.gz && \
	tar -zxvf gsl-$(GSL_VERSION).tar.gz && \
	ln -s gsl-$(GSL_VERSION) gsl

mosek:
	cd external && \
	rm -rf mosek* && \
	wget http://download.mosek.com/stable/$(MOSEK_VERSION)/mosektoolslinux64x86.tar.bz2 && \
	tar jxvf mosektoolslinux64x86.tar.bz2

clean:
	rm -rf build dist *.egg-info

.PHONY: eqtlbma gsl mosek ieqtlbma igsl imosek download install_libs install clean
