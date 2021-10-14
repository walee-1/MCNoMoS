CC=g++

binary1=Maintraj
srcfile=ConfigFiles/config.cpp ConfigFiles/log.cpp
binary2=mainmagJ

EXEC_DIR=bin


all:  $(binary1) $(binary2)

Maintraj: TrajFiles/MaintrajJ.o
	$(shell mkdir -p $(EXEC_DIR))
	$(CC) $^ $(srcfile) -o bin/Maintraj
	cp *.config bin/
	cp *.sh bin/
mainmagJ: MagFiles/mainmagfield3J.o
	$(CC) $^ $(srcfile) -o bin/mainmagJ

clear:
	rm -f TrajFiles/*.o
	rm -f MagFiles/*.o
	rm -f ConfigFiles/*.o

clean: clear
	rm -f $(binary1)
	rm -f $(binary2)
delAll: clean
	$(shell rm -r $(EXEC_DIR))

