RED			=	\033[0;31m

YELLOW		=	\033[0;33m

GREEN		=	\033[0;32m

NC			=	\033[0m

CC			=	g++

NAME		=	energy

CFLAGS		=	-Wall -Wextra -Werror -fopenmp

SRCDIR		=	src

BUILDDIR	=	obj

HDRS		=	headers

SRCS		=	$(wildcard $(SRCDIR)/*.cpp)

OBJS		=	$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SRCS:.cpp=.o))

RM			=	rm -f

all:		$(NAME)

$(NAME):		$(OBJS)
			@echo "$(YELLOW)Linking objects...$(NC)"
			@$(CC) $(CFLAGS) $(OBJS) -o $(NAME)
			@echo "$(GREEN)DONE!$(NC)"

$(OBJS):	|$(BUILDDIR)

$(BUILDDIR):
			@echo "$(YELLOW)Create objects...$(NC)"
			@mkdir $(BUILDDIR)

$(BUILDDIR)/%.o:	$(SRCDIR)/%.cpp $(HDRS)
			$(CC) $(CFLAGS) -I $(HDRS) -c $< -o $@

clean:
			@$(RM) -r $(BUILDDIR)
			@echo "$(RED)Delete objects...$(NC)"

fclean:		clean
			@$(RM) $(NAME)
			@echo "$(RED)Delete executable...$(NC)"
			@echo "$(GREEN)DONE!$(NC)"

re:			fclean all

.PHONY:		all clean fclean re