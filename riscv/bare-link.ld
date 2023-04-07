MEMORY
{
    ram (rwx) : ORIGIN = 0x00000000, LENGTH = 32M /* 32 MB */
}
SECTIONS
{
    . = 0x0;
    .text : {
        *(.text)
    } > ram
    .data : ALIGN(0x1000) {
        *(.data)
    } > ram
    .srodata : {
        *(.srodata)
    } > ram
    .riscv.attributes : {
	*(.riscv.attributes)
    } > ram
    .bss : {
        *(.bss)
    } > ram
}
