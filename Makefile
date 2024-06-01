sync:
	test -d ~/vm/shared && \
		rm -rf ~/vm/shared/voxlap && \
		cp -r . ~/vm/shared/voxlap

.PHONY: sync
