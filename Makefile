voxlap: voxlap5.c sdlmain.c
	gcc $^ \
		-o $@ \
		-m32 \
		-I /home/rbakbashev/build/vox/my/slab6/SDL2inst/include \
		/home/rbakbashev/build/vox/my/slab6/SDL2inst/lib/libSDL2.a \
		-lm \
		-g \
		-funsigned-char
	./$@

clean:
	rm -f voxlap

sync:
	test -d ~/vm/shared && \
		rm -rf ~/vm/shared/voxlap && \
		cp -r . ~/vm/shared/voxlap

.PHONY: clean sync
