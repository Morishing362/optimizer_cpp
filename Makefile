OBJS = $(OUT_DIR)/newton_hvf $(OUT_DIR)/gause_newton_hvf $(OUT_DIR)/LM_hvf $(OUT_DIR)/newton_KandO $(OUT_DIR)/gause_newton_KandO $(OUT_DIR)/LM_KandO $(OUT_DIR)/saikyu_KandO $(OUT_DIR)/test
CC = g++
OUT_DIR = out
SRC_DIR = src

newton_hvf : $(SRC_DIR)/newton_hvf.cpp
	$(CC) -o $(OUT_DIR)/newton_hvf $(SRC_DIR)/newton_hvf.cpp 

gause_newton_hvf : $(SRC_DIR)/gause_newton_hvf.cpp
	$(CC) -o $(OUT_DIR)/gause_newton_hvf $(SRC_DIR)/gause_newton_hvf.cpp 

LM_hvf : $(SRC_DIR)/LM_hvf.cpp
	$(CC) -o $(OUT_DIR)/LM_hvf $(SRC_DIR)/LM_hvf.cpp 

newton_KandO : $(SRC_DIR)/newton_KandO.cpp
	$(CC) -o $(OUT_DIR)/newton_KandO $(SRC_DIR)/newton_KandO.cpp 

gause_newton_KandO : $(SRC_DIR)/gause_newton_KandO.cpp
	$(CC) -o $(OUT_DIR)/gause_newton_KandO $(SRC_DIR)/gause_newton_KandO.cpp 

LM_KandO : $(SRC_DIR)/LM_KandO.cpp
	$(CC) -o $(OUT_DIR)/LM_KandO $(SRC_DIR)/LM_KandO.cpp 

saikyu_KandO : $(SRC_DIR)/saikyu_KandO.cpp
	$(CC) -o $(OUT_DIR)/saikyu_KandO $(SRC_DIR)/saikyu_KandO.cpp 

test : $(SRC_DIR)/test.cpp
	$(CC) -o $(OUT_DIR)/test $(SRC_DIR)/test.cpp 


all :
	make newton_hvf gause_newton_hvf LM_hvf newton_KandO gause_newton_KandO LM_KandO saikyu_KandO test


clean :
	rm -f $(OBJS)