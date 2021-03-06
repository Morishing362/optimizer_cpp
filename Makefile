OBJS = $(OUT_DIR)/newton_hvf \
       $(OUT_DIR)/gauss_newton_hvf \
	   $(OUT_DIR)/LM_hvf \
	   $(OUT_DIR)/newton_KandO \
	   $(OUT_DIR)/gauss_newton_KandO \
	   $(OUT_DIR)/LM_KandO \
	   $(OUT_DIR)/saikyu_KandO \
	   $(OUT_DIR)/newton_Beale \
	   $(OUT_DIR)/gauss_newton_Beale \
	   $(OUT_DIR)/LM_Beale \
	   $(OUT_DIR)/test \
	   $(OUT_DIR)/test_2 \
	   $(OUT_DIR)/Beale_st24 \
	   $(OUT_DIR)/Beale_st25 \
	   $(OUT_DIR)/problem_7_st24 \
	   $(OUT_DIR)/problem_7_st25 \
	   $(OUT_DIR)/LM_problem_7 \
	   $(OUT_DIR)/newton_problem_7 \
	   $(OUT_DIR)/gauss_newton_problem_7
CC = g++
OUT_DIR = out
SRC_DIR = src

newton_hvf : $(SRC_DIR)/newton_hvf.cpp
	$(CC) -o $(OUT_DIR)/newton_hvf $(SRC_DIR)/newton_hvf.cpp 
gauss_newton_hvf : $(SRC_DIR)/gauss_newton_hvf.cpp
	$(CC) -o $(OUT_DIR)/gauss_newton_hvf $(SRC_DIR)/gauss_newton_hvf.cpp 
LM_hvf : $(SRC_DIR)/LM_hvf.cpp
	$(CC) -o $(OUT_DIR)/LM_hvf $(SRC_DIR)/LM_hvf.cpp 


newton_KandO : $(SRC_DIR)/newton_KandO.cpp
	$(CC) -o $(OUT_DIR)/newton_KandO $(SRC_DIR)/newton_KandO.cpp 
gauss_newton_KandO : $(SRC_DIR)/gauss_newton_KandO.cpp
	$(CC) -o $(OUT_DIR)/gauss_newton_KandO $(SRC_DIR)/gauss_newton_KandO.cpp 
LM_KandO : $(SRC_DIR)/LM_KandO.cpp
	$(CC) -o $(OUT_DIR)/LM_KandO $(SRC_DIR)/LM_KandO.cpp 
saikyu_KandO : $(SRC_DIR)/saikyu_KandO.cpp
	$(CC) -o $(OUT_DIR)/saikyu_KandO $(SRC_DIR)/saikyu_KandO.cpp 


newton_Beale : $(SRC_DIR)/newton_Beale.cpp
	$(CC) -o $(OUT_DIR)/newton_Beale $(SRC_DIR)/newton_Beale.cpp
gauss_newton_Beale : $(SRC_DIR)/gauss_newton_Beale.cpp
	$(CC) -o $(OUT_DIR)/gauss_newton_Beale $(SRC_DIR)/gauss_newton_Beale.cpp
LM_Beale : $(SRC_DIR)/LM_Beale.cpp
	$(CC) -o $(OUT_DIR)/LM_Beale $(SRC_DIR)/LM_Beale.cpp
Beale_st24 : $(SRC_DIR)/Beale_st24.cpp
	$(CC) -o $(OUT_DIR)/Beale_st24 $(SRC_DIR)/Beale_st24.cpp
Beale_st25 : $(SRC_DIR)/Beale_st25.cpp
	$(CC) -o $(OUT_DIR)/Beale_st25 $(SRC_DIR)/Beale_st25.cpp

newton_problem_7 : $(SRC_DIR)/newton_problem_7.cpp
	$(CC) -o $(OUT_DIR)/newton_problem_7 $(SRC_DIR)/newton_problem_7.cpp
gauss_newton_problem_7 : $(SRC_DIR)/gauss_newton_problem_7.cpp
	$(CC) -o $(OUT_DIR)/gauss_newton_problem_7 $(SRC_DIR)/gauss_newton_problem_7.cpp
LM_problem_7 : $(SRC_DIR)/LM_problem_7.cpp
	$(CC) -o $(OUT_DIR)/LM_problem_7 $(SRC_DIR)/LM_problem_7.cpp
problem_7_st24 : $(SRC_DIR)/problem_7_st24.cpp
	$(CC) -o $(OUT_DIR)/problem_7_st24 $(SRC_DIR)/problem_7_st24.cpp
problem_7_st25 : $(SRC_DIR)/problem_7_st25.cpp
	$(CC) -o $(OUT_DIR)/problem_7_st25 $(SRC_DIR)/problem_7_st25.cpp


test : $(SRC_DIR)/test.cpp
	$(CC) -o $(OUT_DIR)/test $(SRC_DIR)/test.cpp 
test_2 : $(SRC_DIR)/test_2.cpp
	$(CC) -o $(OUT_DIR)/test_2 $(SRC_DIR)/test_2.cpp


all :
	make $(OBJS)


clean :
	rm -f $(OBJS)