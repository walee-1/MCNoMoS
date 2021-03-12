CC=g++

binaries=MaintrajJ
srcfile=ConfigFiles/config.cpp ConfigFiles/log.cpp

EXEC_DIR=bin


all:  $(binaries)

MaintrajJ: TrajFiles/MaintrajJ.o
	$(shell mkdir -p $(EXEC_DIR))
	$(CC) $^ $(srcfile) -o bin/Maintraj
	cp *.config bin/
	cp *.sh bin/
clear:
	rm -f TrajFiles/*.o
	rm -f MagFiles/*.o
	rm -f ConfigFiles/*.o

clean: clear
	rm -f $(binaries)
delAll: clean
	$(shell rm -r $(EXEC_DIR))

