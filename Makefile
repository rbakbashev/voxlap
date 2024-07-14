CC = gcc

run: voxlap
	./$^

voxlap: voxlap5.c sdlmain.c
	$(CC) \
		$^ \
		-o $@ \
		-lSDL2 \
		-lm \
		-g \
		-Wall \
		-O3 \
		-Wno-misleading-indentation -Wno-char-subscripts \
		-funsigned-char

clean:
	rm -f voxlap

sync:
	test -d ~/vm/shared && \
		rm -rf ~/vm/shared/voxlap && \
		cp -r . ~/vm/shared/voxlap

.PHONY: clean sync
