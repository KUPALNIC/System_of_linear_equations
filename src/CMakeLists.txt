target_sources(${PROJECT_NAME} PRIVATE
    tomas/tomas.cpp
    matrix/matrix.cpp
    matrix/QR/qr.cpp
    matrix/QR/qr_solve.cpp
    CSR/csr.cpp
    vector_oper.cpp
)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)