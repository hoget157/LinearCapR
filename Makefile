PROG := LinearCapR
OBJDIR := .\temp
SRCS := $(wildcard *.cpp) # wildcard関数を用いてファイル内の.cppを全て取得する(配列)
OBJS := $(addprefix $(OBJDIR)\, $(SRCS:%.cpp=%.o)) # %マクロを用いて置換 SRC配列を元に"<ファイル名>.o"の配列を作る。 中間オブジェクトファイル用。
DEPS := $(addprefix $(OBJDIR)\, $(SRCS:%.cpp=%.d)) # %マクロを用いて置換 SRC配列を元に"<ファイル名>.d"の配列を作る。 依存ファイル用。

# 各種設定を変数として定義
CC := g++
CCFLAGS := -O3 -std=c++17 -Wall -pg
# INCLUDEPATH := -I/usr/local/include
# LIBPATH := -L/usr/local/lib
# LIBS := -framework Cocoa -framework OpenGL -lz -ljpeg -lpng

# これが主レシピ
all: $(DEPENDS) $(PROG)

# リンク
$(PROG): $(OBJS)
	$(CC) $(CCFLAGS) -o $@ $^ $(LIBPATH) $(LIBS)

# コンパイル
$(OBJDIR)\\%.o: %.cpp
	$(CC) $(CCFLAGS) $(INCLUDEPATH) -MMD -MP -MF $(<:%.cpp=temp\\%.d) -c $< -o $(<:%.cpp=temp\\%.o)

# "make clean"でターゲットと中間ファイルを消去できるようにする
.PHONY: clean
clean:
	del $(PROG).exe $(OBJS) $(DEPS)

-include $(DEPS)