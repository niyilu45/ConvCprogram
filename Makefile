CC = gcc
OBJECTS = main.o
TARGET = test

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC)  *.c -ggdb3 -pg  -lm  -o  $@

