make _core _cell _geometry _vessel _pde _simulation _mesh _flow _angiogenesis -j 3
export ANGIOGENESIS_MODULE_DIR=./projects/Angiogenesis/src/python/
cp projects/Angiogenesis/_cell.so $ANGIOGENESIS_MODULE_DIR/chaste/population/cell
cp projects/Angiogenesis/_vessel.so $ANGIOGENESIS_MODULE_DIR/chaste/population/vessel
cp projects/Angiogenesis/_mesh.so $ANGIOGENESIS_MODULE_DIR/chaste/mesh
cp projects/Angiogenesis/_pde.so $ANGIOGENESIS_MODULE_DIR/chaste/pde
cp projects/Angiogenesis/_simulation.so $ANGIOGENESIS_MODULE_DIR/chaste/simulation
cp projects/Angiogenesis/_flow.so $ANGIOGENESIS_MODULE_DIR/chaste/simulation
cp projects/Angiogenesis/_angiogenesis.so $ANGIOGENESIS_MODULE_DIR/chaste/simulation
cp projects/Angiogenesis/_core.so $ANGIOGENESIS_MODULE_DIR/chaste/core
cp projects/Angiogenesis/_geometry.so $ANGIOGENESIS_MODULE_DIR/chaste/geometry