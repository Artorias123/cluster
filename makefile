CC	=g++
CFLAGS	= -Wall -g
INCLUDE	= ./include
OBJPATH = obj
SRCPATH = src
TARGET	= run
OBJS    = $(addprefix obj/,$(patsubst %.cpp,%.o,$(notdir $(SRCS))))
SRCS    = $(wildcard src/*.cpp)
vpath %.cpp src
vpath %.o obj
$(TARGET):		$(OBJS)
			$(CC) -o $@ $^ -I$(INCLUDE)
$(OBJPATH)/%.o:		$(SRCPATH)/%.cpp
			$(CC) $(CFLAGS) -c $^ -I$(INCLUDE) -o $@
clean:			
			rm obj/*
			rm test
