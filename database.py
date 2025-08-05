import os
import typing
from typing import Literal, overload

import numpy as np
import tables as tbl
from scipy.fft import rfft

from util import assay_index_to_type

COMPRESSION_LIB = "blosc2:zstd"
COMPRESSION_LEVEL = 1
COMPRESSION_FILTER = tbl.Filters(
    complib=COMPRESSION_LIB,  # pyright: ignore [reportArgumentType]
    complevel=COMPRESSION_LEVEL,
    shuffle=True,
)
NUM_EXPECTED_ROWS = 3_000
MAX_RAW_SAMPLES = 2702
NUM_SAMPLES = 2400
NUM_SAMPLES_FFT = 1201
MAX_ORDER = 2


class ParticipantIDMapDescription(tbl.IsDescription):
    participant_id = tbl.StringCol(16)  # pyright: ignore [reportCallIssue]


class CPHS2017TrialMetadataDescription(tbl.IsDescription):
    timestamp = tbl.Time64Col()  # pyright: ignore [reportCallIssue]
    protocol_name = tbl.StringCol(16)  # pyright: ignore [reportCallIssue]
    participant_num = tbl.UInt8Col()  # pyright: ignore [reportCallIssue]
    participant_id = tbl.StringCol(16)  # pyright: ignore [reportCallIssue]
    order = tbl.UInt8Col()  # pyright: ignore [reportCallIssue]
    num_reps = tbl.UInt8Col()  # pyright: ignore [reportCallIssue]
    reference_type = tbl.Int8Col()  # pyright: ignore [reportCallIssue]
    disturbance_type = tbl.Int8Col()  # pyright: ignore [reportCallIssue]
    random_seed = tbl.UInt32Col()  # pyright: ignore [reportCallIssue]
    phase_shift_id = tbl.UInt8Col()  # pyright: ignore [reportCallIssue]
    assay_type = tbl.UInt8Col()  # pyright: ignore [reportCallIssue]
    assay_index = tbl.UInt8Col()  # pyright: ignore [reportCallIssue]
    sign = tbl.Int8Col()  # pyright: ignore [reportCallIssue]


CPHS2017DataKey = Literal[
    "timestamp",
    "protocol_name",
    "participant_id",
    "participant_num",
    "order",
    "num_reps",
    "reference_type",
    "disturbance_type",
    "random_seed",
    "phase_shift_id",
    "assay_type",
    "assay_index",
    "sign",
]


class CPHS2017DatabaseRow:
    def __init__(self, row):
        self._row = row

    def __len__(self) -> int:
        return len(self._row)  # pyright: ignore [reportArgumentType]

    @property
    def data(self):
        return self._row

    @overload
    def __getitem__(
        self, key: Literal["timestamp"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.float64]]:
        return self._row[key]

    @overload
    def __getitem__(
        self, key: Literal["protocol_name"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.str_]]:
        return self._row[key]

    @overload
    def __getitem__(
        self, key: Literal["participant_id"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.str_]]:
        return self._row[key]

    @overload
    def __getitem__(
        self, key: Literal["participant_num"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row[key]

    @overload
    def __getitem__(
        self, key: Literal["order"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row[key]

    @overload
    def __getitem__(
        self, key: Literal["num_reps"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row[key]

    @overload
    def __getitem__(
        self, key: Literal["reference_type"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row[key]

    @overload
    def __getitem__(
        self, key: Literal["disturbance_type"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row[key]

    @overload
    def __getitem__(
        self, key: Literal["random_seed"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.uint32]]:
        return self._row[key]

    @overload
    def __getitem__(
        self, key: Literal["phase_shift_id"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row[key]

    @overload
    def __getitem__(
        self, key: Literal["assay_type"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row[key]

    @overload
    def __getitem__(
        self, key: Literal["assay_index"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row[key]

    @overload
    def __getitem__(
        self, key: Literal["sign"]
    ) -> np.ndarray[tuple[int, ...], np.dtype[np.int8]]:
        return self._row[key]

    def __getitem__(
        self, key: CPHS2017DataKey
    ) -> (
        np.ndarray[tuple[int, ...], np.dtype[np.uint8]]
        | np.ndarray[tuple[int, ...], np.dtype[np.int8]]
        | np.ndarray[tuple[int, ...], np.dtype[np.uint32]]
        | np.ndarray[tuple[int, ...], np.dtype[np.float64]]
        | np.ndarray[tuple[int, ...], np.dtype[np.str_]]
    ): ...
    @property
    def timestamp(self) -> np.ndarray[tuple[int, ...], np.dtype[np.float64]]:
        return self._row["timestamp"]

    @timestamp.setter
    def timestamp(self, value: float) -> None:
        self._row["timestamp"] = value  # Implicitly casts float → np.float64

    @property
    def protocol_name(self) -> np.ndarray[tuple[int, ...], np.dtype[np.str_]]:
        return self._row["protocol_name"]

    @protocol_name.setter
    def protocol_name(self, value: str) -> None:
        self._row["protocol_name"] = value  # Implicitly casts str → np.str_

    @property
    def participant_id(self) -> np.ndarray[tuple[int, ...], np.dtype[np.str_]]:
        return self._row["participant_id"]

    @participant_id.setter
    def participant_id(self, value: str) -> None:
        self._row["participant_id"] = value  # Implicitly casts str → np.str_

    @property
    def participant_num(self) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row["participant_num"]

    @participant_num.setter
    def participant_num(self, value: int) -> None:
        self._row["participant_num"] = value  # Implicitly casts int → np.uint8

    @property
    def order(self) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row["order"]

    @order.setter
    def order(self, value: int) -> None:
        self._row["order"] = value  # Implicitly casts int → np.uint8

    @property
    def num_reps(self) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row["num_reps"]

    @num_reps.setter
    def num_reps(self, value: int) -> None:
        self._row["num_reps"] = value  # Implicitly casts int → np.uint8

    @property
    def reference_type(self) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row["reference_type"]

    @reference_type.setter
    def reference_type(self, value: int) -> None:
        self._row["reference_type"] = value  # Implicitly casts int → np.uint8

    @property
    def disturbance_type(self) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row["disturbance_type"]

    @disturbance_type.setter
    def disturbance_type(self, value: int) -> None:
        self._row["disturbance_type"] = value  # Implicitly casts int → np.uint8

    @property
    def random_seed(self) -> np.ndarray[tuple[int, ...], np.dtype[np.uint32]]:
        return self._row["random_seed"]

    @random_seed.setter
    def random_seed(self, value: int) -> None:
        self._row["random_seed"] = value  # Implicitly casts int → np.uint32

    @property
    def phase_shift_id(self) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row["phase_shift_id"]

    @phase_shift_id.setter
    def phase_shift_id(self, value: int) -> None:
        self._row["phase_shift_id"] = value  # Implicitly casts int → np.uint8

    @property
    def assay_type(self) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row["assay_type"]

    @assay_type.setter
    def assay_type(self, value: int) -> None:
        self._row["assay_type"] = value  # Implicitly casts int → np.uint8

    @property
    def assay_index(self) -> np.ndarray[tuple[int, ...], np.dtype[np.uint8]]:
        return self._row["assay_index"]

    @assay_index.setter
    def assay_index(self, value: int) -> None:
        self._row["assay_index"] = value  # Implicitly casts int → np.uint8

    @property
    def sign(self) -> np.ndarray[tuple[int, ...], np.dtype[np.int8]]:
        return self._row["sign"]

    @sign.setter
    def sign(self, value: int) -> None:
        self._row["sign"] = value  # Implicitly casts int → np.uint8


class CPHS2017Database:
    def __init__(self, name: str = "cphs2017_pytable.h5") -> None:
        self._database_file = self._ensure_initialized(name)
        self._metadata = typing.cast(
            tbl.Table, self._database_file.get_node("/metadata")
        )

        self._participant_id_tbl = typing.cast(
            tbl.Table, self._database_file.get_node("/participant_id_table")
        )

        self._participant_id_map = {
            id: i for i, id in self._participant_id_tbl.col("participant_id")
        }

        self._raw_time = typing.cast(
            tbl.EArray, self._database_file.get_node("/raw_time")
        )
        self._raw_real_time = typing.cast(
            tbl.EArray, self._database_file.get_node("/raw_real_time")
        )
        self._raw_reference = typing.cast(
            tbl.EArray, self._database_file.get_node("/raw_reference")
        )
        self._raw_disturbance = typing.cast(
            tbl.EArray, self._database_file.get_node("/raw_disturbance")
        )
        self._raw_input = typing.cast(
            tbl.EArray, self._database_file.get_node("/raw_input")
        )
        self._raw_state = typing.cast(
            tbl.EArray, self._database_file.get_node("/raw_state")
        )
        self._time = typing.cast(tbl.EArray, self._database_file.get_node("/time"))
        self._real_time = typing.cast(
            tbl.EArray, self._database_file.get_node("/real_time")
        )
        self._reference = typing.cast(
            tbl.EArray, self._database_file.get_node("/reference")
        )
        self._disturbance = typing.cast(
            tbl.EArray, self._database_file.get_node("/disturbance")
        )
        self._input = typing.cast(tbl.EArray, self._database_file.get_node("/input"))
        self._state = typing.cast(tbl.EArray, self._database_file.get_node("/state"))

        self._fft_reference = typing.cast(
            tbl.EArray, self._database_file.get_node("/fft_reference")
        )
        self._fft_disturbance = typing.cast(
            tbl.EArray, self._database_file.get_node("/fft_disturbance")
        )
        self._fft_input = typing.cast(
            tbl.EArray, self._database_file.get_node("/fft_input")
        )

        self._tf_disturbance_to_input = typing.cast(
            tbl.EArray, self._database_file.get_node("/tf_disturbance_to_input")
        )
        self._tf_reference_to_input = typing.cast(
            tbl.EArray, self._database_file.get_node("/tf_reference_to_input")
        )

    @property
    def raw_time(self) -> tbl.EArray:
        """Raw time data array (read-only)"""
        return self._raw_time

    @property
    def raw_real_time(self) -> tbl.EArray:
        """Raw real time data array (read-only)"""
        return self._raw_real_time

    @property
    def raw_reference(self) -> tbl.EArray:
        """Raw reference signal data array (read-only)"""
        return self._raw_reference

    @property
    def raw_disturbance(self) -> tbl.EArray:
        """Raw disturbance signal data array (read-only)"""
        return self._raw_disturbance

    @property
    def raw_input(self) -> tbl.EArray:
        """Raw input signal data array (read-only)"""
        return self._raw_input

    @property
    def raw_state(self) -> tbl.EArray:
        """Raw state data array (read-only)"""
        return self._raw_state

    @property
    def time(self) -> tbl.EArray:
        """Processed time data array (read-only)"""
        return self._time

    @property
    def real_time(self) -> tbl.EArray:
        """Processed real time data array (read-only)"""
        return self._real_time

    @property
    def reference(self) -> tbl.EArray:
        """Processed reference signal data array (read-only)"""
        return self._reference

    @property
    def disturbance(self) -> tbl.EArray:
        """Processed disturbance signal data array (read-only)"""
        return self._disturbance

    @property
    def input(self) -> tbl.EArray:
        """Processed input signal data array (read-only)"""
        return self._input

    @property
    def state(self) -> tbl.EArray:
        """Processed state data array (read-only)"""
        return self._state

    @property
    def fft_reference(self) -> tbl.EArray:
        """FFT of the reference signal (read-only)"""
        return self._fft_reference

    @property
    def fft_disturbance(self) -> tbl.EArray:
        """FFT of the disturbance signal (read-only)"""
        return self._fft_disturbance

    @property
    def fft_input(self) -> tbl.EArray:
        """FFT of the input signal (read-only)"""
        return self._fft_input

    @property
    def tf_disturbance_to_input(self) -> tbl.EArray:
        """Transfer function from disturbance to input (read-only)"""
        return self._tf_disturbance_to_input

    @property
    def tf_reference_to_input(self) -> tbl.EArray:
        """Transfer function from reference to input (read-only)"""
        return self._tf_reference_to_input

    def append(
        self,
        *,
        timestamp: float,
        protocol_name: str,
        participant_id: str,
        order: int,
        num_reps: int,
        reference_type: int,
        disturbance_type: int,
        random_seed: int,
        phase_shift_id: int,
        assay_index: int,
        sign: int,
        raw_time: list[float] | np.ndarray[tuple[int, ...], np.dtype[np.float64]],
        raw_real_time: list[float] | np.ndarray[tuple[int, ...], np.dtype[np.float64]],
        raw_input: list[float] | np.ndarray[tuple[int, ...], np.dtype[np.float64]],
        raw_reference: list[float] | np.ndarray[tuple[int, ...], np.dtype[np.float64]],
        raw_disturbance: list[float]
        | np.ndarray[tuple[int, ...], np.dtype[np.float64]],
        raw_state: list[list[float]]
        | np.ndarray[tuple[int, ...], np.dtype[np.float64]],
    ) -> None:
        row = CPHS2017DatabaseRow(self._metadata.row)

        row.timestamp = timestamp
        row.protocol_name = protocol_name
        row.participant_id = participant_id

        if participant_id not in self._participant_id_map:
            self._participant_id_map[participant_id] = self._participant_id_tbl.nrows
            participant_id_row = self._participant_id_tbl.row
            participant_id_row['participant_id'] = participant_id

            participant_id_row.append()
            self._participant_id_tbl.flush()

        row.participant_num = self._participant_id_map[participant_id]
        row.order = order
        row.num_reps = num_reps
        row.reference_type = reference_type
        row.disturbance_type = disturbance_type
        row.random_seed = random_seed
        row.phase_shift_id = phase_shift_id
        row.assay_type = assay_index_to_type(assay_index)
        row.assay_index = assay_index
        row.sign = sign
        row.data.append()

        raw_time_arr = np.asarray(raw_time)
        N_raw = raw_time_arr.shape[0]
        self._raw_time.append(raw_time_arr[np.newaxis, :])

        raw_real_time_arr = np.asarray(raw_real_time)
        self._raw_real_time.append(raw_real_time_arr[np.newaxis, :])

        raw_reference_arr = np.asarray(raw_reference)
        self._raw_reference.append(raw_reference_arr[np.newaxis, :])

        raw_disturbance_arr = np.asarray(raw_disturbance)
        self._raw_disturbance.append(raw_disturbance_arr[np.newaxis, :])

        raw_input_arr = np.asarray(raw_input)
        self._raw_input.append(raw_input_arr[np.newaxis, :])

        raw_state_arr = np.asarray(raw_state).reshape((MAX_RAW_SAMPLES, -1))
        x_raw = np.pad(
            raw_state_arr,
            pad_width=(
                (0, MAX_RAW_SAMPLES - N_raw),
                (0, MAX_ORDER - raw_state_arr.shape[1]),
            ),
            mode="constant",
            constant_values=np.nan,
        )
        self._raw_state.append(x_raw[np.newaxis, :])

        self._time.append(raw_time_arr[np.newaxis, -NUM_SAMPLES:])
        self._real_time.append(raw_real_time_arr[np.newaxis, -NUM_SAMPLES:])

        reference_arr = raw_reference_arr[-NUM_SAMPLES:]
        self._reference.append(reference_arr[np.newaxis, :])

        disturbance_arr = raw_disturbance_arr[-NUM_SAMPLES:]
        self._disturbance.append(disturbance_arr[np.newaxis, :])

        input_arr = raw_input_arr[-NUM_SAMPLES:]
        self._input.append(input_arr[np.newaxis, :])

        self._state.append(x_raw[np.newaxis, -NUM_SAMPLES:])

        R = rfft(reference_arr)[np.newaxis, :]
        self._fft_reference.append(R)

        D = rfft(disturbance_arr)[np.newaxis, :]
        self._fft_disturbance.append(D)

        U = rfft(input_arr)[np.newaxis, :]
        self._fft_input.append(U)

        self._tf_reference_to_input.append(U / R)
        self._tf_disturbance_to_input.append(U / D)

    def __iter__(self):
        return map(CPHS2017DatabaseRow, iter(self._metadata))

    def __enter__(self, name: str = "data.h5"):
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.close()

    def _create_database(self, f: tbl.File):
        f.create_table(
            "/",
            "metadata",
            CPHS2017TrialMetadataDescription,  # pyright: ignore [reportArgumentType]
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_table(
            "/",
            "participant_id_table",
            ParticipantIDMapDescription,
            expectedrows=20,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "raw_time",
            tbl.Float64Atom(),
            shape=(
                0,
                MAX_RAW_SAMPLES,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "raw_real_time",
            tbl.Float64Atom(),
            shape=(
                0,
                MAX_RAW_SAMPLES,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "raw_reference",
            tbl.Float64Atom(),
            shape=(
                0,
                MAX_RAW_SAMPLES,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "raw_disturbance",
            tbl.Float64Atom(),
            shape=(
                0,
                MAX_RAW_SAMPLES,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "raw_input",
            tbl.Float64Atom(),
            shape=(
                0,
                MAX_RAW_SAMPLES,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "raw_state",
            tbl.Float64Atom(),
            shape=(0, MAX_RAW_SAMPLES, MAX_ORDER),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "time",
            tbl.Float64Atom(),
            shape=(
                0,
                NUM_SAMPLES,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "real_time",
            tbl.Float64Atom(),
            shape=(
                0,
                NUM_SAMPLES,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "reference",
            tbl.Float64Atom(),
            shape=(
                0,
                NUM_SAMPLES,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "disturbance",
            tbl.Float64Atom(),
            shape=(
                0,
                NUM_SAMPLES,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "input",
            tbl.Float64Atom(),
            shape=(
                0,
                NUM_SAMPLES,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "state",
            tbl.Float64Atom(),
            shape=(0, NUM_SAMPLES, MAX_ORDER),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "fft_input",
            tbl.ComplexAtom(itemsize=16),
            shape=(
                0,
                NUM_SAMPLES_FFT,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "fft_reference",
            tbl.ComplexAtom(itemsize=16),
            shape=(
                0,
                NUM_SAMPLES_FFT,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "fft_disturbance",
            tbl.ComplexAtom(itemsize=16),
            shape=(
                0,
                NUM_SAMPLES_FFT,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "tf_disturbance_to_input",
            tbl.ComplexAtom(itemsize=16),
            shape=(
                0,
                NUM_SAMPLES_FFT,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

        f.create_earray(
            "/",
            "tf_reference_to_input",
            tbl.ComplexAtom(itemsize=16),
            shape=(
                0,
                NUM_SAMPLES_FFT,
            ),
            expectedrows=NUM_EXPECTED_ROWS,
            filters=COMPRESSION_FILTER,
        )

    def _ensure_initialized(self, name: str = "cphs2017.h5") -> tbl.File:
        if os.path.exists(name):
            print(f"Found database {name}")

            f = tbl.open_file(name, "a", title="Test file")

            if "/metadata" not in f:
                print(f"Logs table not found in database {name}. Creating...")

                self._create_database(f)

                print(f"Created logs table in database {name}.")
                return f

            print(f"Logs table found in database {name}.")
            return f

        print(f"Database {name} not found. Creating...")
        f = tbl.open_file(name, "a", title="Test file")
        self._create_database(f)

        print(f"Created {name}.")

        return f

    def close(self):
        self._metadata.close()
        self._participant_id_tbl.close()

        self._raw_reference.close()
        self._raw_disturbance.close()
        self._raw_input.close()
        self._raw_state.close()

        self._reference.close()
        self._disturbance.close()
        self._input.close()
        self._state.close()

        self._fft_reference.close()
        self._fft_disturbance.close()
        self._fft_input.close()

    def get_where_list(
        self,
        condition: str,
        condvars: dict[str, tbl.Column | np.ndarray] | None = None,
        sort: bool = False,
        start: str | None = None,
        stop: str | None = None,
        step: str | None = None,
    ) -> list[int]:
        r"""Iterate over values fulfilling a condition.

        This method returns a Row iterator (see
        :ref:`RowClassDescr`) which only selects rows in the
        table that satisfy the given condition (an expression-
        like string).

        The condvars mapping may be used to define the variable
        names appearing in the condition. condvars should consist
        of identifier-like strings pointing to Column (see
        :ref:`ColumnClassDescr`) instances *of this table*, or to
        other values (which will be converted to arrays). A
        default set of condition variables is provided where each
        top-level, non-nested column with an identifier-like name
        appears. Variables in condvars override the default ones.

        When condvars is not provided or None, the current local
        and global namespace is sought instead of condvars. The
        previous mechanism is mostly intended for interactive
        usage. To disable it, just specify a (maybe empty)
        mapping as condvars.

        If a range is supplied (by setting some of the start,
        stop or step parameters), only the rows in that range and
        fulfilling the condition are used. The meaning of the
        start, stop and step parameters is the same as for Python
        slices.

        When possible, indexed columns participating in the
        condition will be used to speed up the search. It is
        recommended that you place the indexed columns as left
        and out in the condition as possible. Anyway, this method
        has always better performance than regular Python
        selections on the table.

        You can mix this method with regular Python selections in
        order to support even more complex queries. It is
        strongly recommended that you pass the most restrictive
        condition as the parameter to this method if you want to
        achieve maximum performance.

        .. warning::

        When in the middle of a table row iterator, you should
        not use methods that can change the number of rows in the
        table (like :meth:`Table.append` or
        :meth:`Table.remove_rows`) or unexpected errors will
        happen.

        Examples --------

        ::

        passvalues = [     row["col3"]     for row in
        table.where(         "(col1 > 0) & (col2 <= 20)",
        step=COMPRESSION_LEVEL )     if
        your_function(row["col2"]) ] print("Values that pass the
        cuts:", passvalues)

        .. note::

        A special care should be taken when the query condition
        includes string literals.

        Let's assume that the table ``table`` has the following
        structure::

        class Record(IsDescription):     col1 = StringCol( 4 )  #
        4-character String of bytes     col2 = IntCol() col3 =
        FloatCol()

        The type of "col1" corresponds to strings of bytes.

        Any condition involving "col1" should be written using
        the appropriate type for string literals in order to
        avoid :exc:`TypeError`\ s.

        The code below will fail with a :exc:`TypeError`::

        condition = 'col1 == "AAAA"' for record in
        table.where(condition):  # TypeError in Python3     # do
        something with "record"

        The reason is that in Python 3 "condition" implies a
        comparison between a string of bytes ("col1" contents)
        and a unicode literal ("AAAA").

        The correct way to write the condition is::

        condition = 'col1 == b"AAAA"'

        .. versionchanged:: 3.0    The start, stop and step
        parameters now behave like in slice.

        """

        return typing.cast(
            list[int],
            self._metadata.get_where_list(condition, condvars, sort, start, stop, step),
        )
