#!/bin/bash

is_readable() {

    if [[ ! -f "$1" ]]; then
        echo "$1 not found."
        return
    else
        echo >> $1
        sed -l'a\' $1
    fi

    if [[ ! -r "$1" ]]; then
        echo "$1 is not readable."
        if [ -n "$(tail -c 1 "$1")" ]; then
            echo "" >> $1
            echo "new_line appended"
        fi
        if [[  -r "$1" ]]; then
            echo "$1 is readable."
        else
            sed -i '$a\' $1
            echo "new_line sed"
        fi
    fi
}

INPUT_FILE="$1"
OUTPUT_FILE="$1.new"
CHUNK_SIZE=100

# Check if file exists

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "$INPUT_FILE not found."
    exit 0
fi

is_readable $INPUT_FILE
command_failed=0

if ! line_count=$(grep -c '$' "$INPUT_FILE"); then
    echo "Grep failed!"
    command_failed=1
fi

if [ ${command_failed:-0} -eq 1 ]; then
    fold -w "$CHUNK_SIZE" "$INPUT_FILE" > "$OUTPUT_FILE" || true
    is_readable $OUTPUT_FILE
fi

if [[  -f "$OUTPUT_FILE" ]]; then
    if [[  -r "$OUTPUT_FILE" ]]; then
        mv $OUTPUT_FILE $INPUT_FILE
        touch "$INPUT_FILE".solved
        exit 0
    fi
fi

echo "Process incomplete. check $INPUT_FILE $OUTPUT_FILE"