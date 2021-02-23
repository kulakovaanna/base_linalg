cmake_minimum_required(VERSION 2.8)

project(base_linalg LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCE_EXE main.cpp)		              # Установка переменной со списком исходников для исполняемого файла
set(SOURCE_LIB linalg.cpp)			          # Тоже самое, но для библиотеки
add_library(linalg STATIC ${SOURCE_LIB})	# Создание статической библиотеки

add_executable(main ${SOURCE_EXE})	      # Создает исполняемый файл с именем main
target_link_libraries(main linalg)		    # Линковка программы с библиотекой
